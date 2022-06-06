__author__ = "Yiyang Fu <yiyang.fu@mail.utoronto.ca>"
__date__ = "2022-04-14"

from ..NSfracStep import *
from ..problem_base import *
import numpy as np
from os import getcwd, makedirs
import pickle
import mpi4py
from math import exp, cos, sin, pi, sqrt

# Override some problem specific parameters
def problem_parameters(commandline_kwargs, NS_parameters, NS_expressions, **NS_namespace):
    if "restart_folder" in commandline_kwargs.keys():
        restart_folder = commandline_kwargs["restart_folder"]
        restart_folder = path.join(getcwd(), restart_folder)
        f = open(path.join(path.dirname(path.abspath(__file__)), restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        f.close()
        NS_parameters['restart_folder'] = restart_folder
        globals().update(NS_parameters)
    else:
        NS_parameters.update(
            nu=0.037736,
            T=0.1,
            dt=0.0005,
            folder="folder_name",
            velocity_degree=2,
            pressure_degree=1,
            plot_interval=100000,
            save_step=1,
            checkpoint=10,
            print_intermediate_info=100000,
            use_krylov_solvers=True)
        NS_parameters['krylov_solvers'] = {'monitor_convergence': False,
                                          'report': False,
                                          'error_on_nonconvergence': True,
                                          'nonzero_initial_guess':True,
                                          'maximum_iterations':300,
                                          'relative_tolerance': 1e-5,
                                          'absolute_tolerance': 1e-5}
    
    # NS_parameters['velocity_krylov_solver']={
    #                                         'solver_type':'gmres',
    #                                         # 'solver_type':'bicgstab',
    #                                         'preconditioner_type':'hypre_amg'}
    #                                         # 'preconditioner_type':'ilu'}

def pre_solve_hook(restart_folder, **NS_namespace):
    if restart_folder is None:
        d = dict(Q_bca_1=Q_bca_1,Q_lcc_1=Q_lcc_1,Q_lsub_1=Q_lsub_1,Q_aao_1=Q_aao_1,rcr_aao=rcr_aao,rcr_bca=rcr_bca,rcr_lsub=rcr_lsub,rcr_lcc=rcr_lcc)
        return d
    else:
        return {}

# Specify boundary conditions
def create_bcs(V,Q, **NS_namespace):
    flow_rate = 1.0 #unit: cc
    bcp_outflow = DirichletBC(Q, Constant(0), capid_markers, capid_seq['outlet'])
    bcp_outflow_BCA = DirichletBC(Q, Constant(0), capid_markers, capid_seq['BCA'])
    bcp_outflow_LCC = DirichletBC(Q, Constant(0), capid_markers, capid_seq['LCC'])
    bcp_outflow_LSUB = DirichletBC(Q, Constant(0), capid_markers, capid_seq['LSUB'])
    ds = Measure('ds', domain=mesh, subdomain_data=capid_markers)
    inlet_area=assemble(Constant(1)*ds(capid_seq['inlet']))
    print('Inlet area of aorta: ', inlet_area)
    insc_R = 1.77672
    inscrib_circ_area = pi*insc_R*insc_R
    v_max =(2*flow_rate/(inscrib_circ_area))
 
    nx = -0.49872
    ny = -0.226005
    nz = 0.836779
    inlet_centre = [0.0912544, 1.68265, 182.677]
    cx = inlet_centre[0]
    cy = inlet_centre[1]
    cz = inlet_centre[2]
    dis_sk_c = sqrt((sk_x-cx)*(sk_x-cx)+(sk_y-cy)*(sk_y-cy)+(sk_z-cz)*(sk_z-cz))
    u_in_pois_x = Expression(("(nx*v)*(1.0 - (pow(x[0]-cx, 2) + pow(x[1]-cy, 2) + pow(x[2] - cz, 2))/pow(alpha, 2)"), cx=cx, cy=cy, cz=cz, nx=nx, v=v_max, alpha=insc_R,degree=2)
    u_in_pois_y = Expression(("(ny*v)*(1.0 - (pow(x[0]-cx, 2) + pow(x[1]-cy, 2) + pow(x[2] - cz, 2))/pow(alpha, 2)"), cx=cx, cy=cy, cz=cz, ny=ny, v=v_max, alpha=insc_R,degree=2)
    u_in_pois_z = Expression(("(nz*v)*(1.0 - (pow(x[0]-cx, 2) + pow(x[1]-cy, 2) + pow(x[2] - cz, 2))/pow(alpha, 2)"), cx=cx, cy=cy, cz=cz, nz=nz, v=v_max, alpha=insc_R,degree=2)

    bcu_inflow_x = DirichletBC(V, u_in_pois_x, capid_markers, capid_seq['inlet'])
    bcu_inflow_y = DirichletBC(V, u_in_pois_y, capid_markers, capid_seq['inlet'])
    bcu_inflow_z = DirichletBC(V, u_in_pois_z, capid_markers, capid_seq['inlet'])
    bcu_walls = DirichletBC(V, Constant(0), capid_markers, 1)
    return dict(u0=[bcu_inflow_x, bcu_walls],
                u1=[bcu_inflow_y, bcu_walls],
                u2=[bcu_inflow_z, bcu_walls],
                p=[bcp_outflow, bcp_outflow_BCA, bcp_outflow_LCC, bcp_outflow_LSUB])


def initialize(x_, x_1, x_2, bcs,restart_folder,newfolder, tstepfiles, tstep, Q_bca_1,Q_lcc_1,Q_lsub_1,Q_aao_1,rcr_aao,rcr_bca,rcr_lsub,rcr_lcc,**NS_namespace):
    add_function_to_tstepfiles(shear_stress_, newfolder, tstepfiles, tstep)
    if restart_folder is None:
        for ui in x_1:
            [bc.apply(x_1[ui]) for bc in bcs[ui]]
        for ui in x_2:
            [bc.apply(x_2[ui]) for bc in bcs[ui]]
    else:
        comm = mpi4py.MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank==0:
            print('---Restarted variables from last timestep---')
            print('current hook: initialize')
            print('Q_bca_1:',Q_bca_1)
            print('Q_lcc_1:',Q_lcc_1)
            print('Q_lsub_1:',Q_lsub_1)
            print('Q_aao_1:',Q_aao_1)
            print('rcr_aao:',rcr_aao)
            print('rcr_bca:',rcr_bca)
            print('rcr_lcc:',rcr_lcc)
            print('rcr_lsub:',rcr_lsub)


#backflow stabilization for high-performance IPCS-ABCN solver
def velocity_tentative_hook(ui,A,u_ab,u,v,q_1,b,x_1, **NS_namespace):

    ds = Measure('ds', domain=mesh, subdomain_data=capid_markers)
    n_normal = FacetNormal(mesh)
    beta = 0.2
    #matrix for stabilization term
    ST1 = assemble(dotun_(u_ab,n_normal)*dot(u, v)*ds(capid_seq['outlet']))
    # add stablization term to variational form for aorta outlet
    # print("Adding stablization Term to bilinear form matrix A")
    A.axpy(-beta/2,ST1, True)
    # add stabilization term to rhs linear form
    b[ui].axpy(beta/2, ST1 * x_1[ui])
    # for BCA
    ST2 = assemble(dotun_(u_ab,n_normal)*dot(u, v)*ds(capid_seq['BCA']))
    A.axpy(-beta/2,ST2, True)
    b[ui].axpy(beta/2, ST2 * x_1[ui])
    # for LSUB
    ST3 = assemble(dotun_(u_ab,n_normal)*dot(u, v)*ds(capid_seq['LSUB']))
    A.axpy(-beta/2,ST3, True)
    b[ui].axpy(beta/2, ST3 * x_1[ui])
    # for LCC
    ST4 = assemble(dotun_(u_ab,n_normal)*dot(u, v)*ds(capid_seq['LCC']))
    A.axpy(-beta/2,ST4, True)
    b[ui].axpy(beta/2, ST4 * x_1[ui])

def temporal_hook(t, dt, x_, u_, p_, tstep, save_step, V, Q, bcs,Q_bca_1,Q_lcc_1,Q_lsub_1,Q_aao_1,rcr_aao,rcr_bca,rcr_lsub,rcr_lcc, NS_parameters, **NS_namespace):
    # if tstep % plot_interval == 0 and not testing:
    ds = Measure('ds', domain=mesh, subdomain_data=capid_markers)
    n_normal = FacetNormal(mesh)
    outlet_area=assemble(Constant(1)*ds(capid_seq['outlet']))
    flux = dot(u_, n_normal)*ds(capid_seq['outlet'])
    Q_flux = assemble(flux)
    flux_inlet = dot(u_, n_normal)*ds(capid_seq['inlet'])
    Q_inlet = assemble(flux_inlet)
    # calcualte pressure at aorta outlet
    p_value = Q_flux*Rp_aao
    # update the rcr term
    temp1 = 1/2*dt*(exp(-dt/(Rd_aao*Cap_aao))/Cap_aao*Q_aao_1[0] + Q_flux/Cap_aao)
    temp2 = exp(-dt/(Rd_aao*Cap_aao))*rcr_aao[0]
    temp = temp1+temp2
    p_value += temp
    rcr_aao[0] = temp
    Q_aao_1[0] = Q_flux
    #update bcp at aorta outlet
    bcp_outflow = DirichletBC(Q, Constant(p_value), capid_markers, capid_seq['outlet'])
    #calculate outlets area and flux
    outlet_area_BCA=assemble(Constant(1)*ds(capid_seq['BCA']))
    outlet_area_LSUB=assemble(Constant(1)*ds(capid_seq['LSUB']))
    outlet_area_LCC=assemble(Constant(1)*ds(capid_seq['LCC']))
    flux_BCA = dot(u_, n_normal)*ds(capid_seq['BCA'])
    flux_LSUB = dot(u_, n_normal)*ds(capid_seq['LSUB'])
    flux_LCC = dot(u_, n_normal)*ds(capid_seq['LCC'])
    Q_flux_BCA = assemble(flux_BCA)
    Q_flux_LSUB = assemble(flux_LSUB)
    Q_flux_LCC = assemble(flux_LCC)
    # calculate p for BCA, add rcr term to pressure at BCA
    p_value_BCA = Q_flux_BCA*Rp_BCA
    temp1_bca = 1/2*dt*(exp(-dt/(Rd_BCA*Cap_BCA))/Cap_BCA*Q_bca_1[0] + Q_flux_BCA/Cap_BCA)
    temp2_bca = exp(-dt/(Rd_BCA*Cap_BCA))*rcr_bca[0]
    temp_bca = temp1_bca + temp2_bca
    p_value_BCA += temp_bca
    rcr_bca[0] = temp_bca
    Q_bca_1[0] = Q_flux_BCA
    # calculate p for LSUB, add rcr term to pressure at LSUB
    p_value_LSUB = Q_flux_LSUB*Rp_LSUB
    temp1_lsub = 1/2*dt*(exp(-dt/(Rd_LSUB*Cap_LSUB))/Cap_LSUB*Q_lsub_1[0] + Q_flux_LSUB/Cap_LSUB)
    temp2_lsub = exp(-dt/(Rd_LSUB*Cap_LSUB))*rcr_lsub[0]
    temp_lsub = temp1_lsub + temp2_lsub
    p_value_LSUB += temp_lsub
    rcr_lsub[0] = temp_lsub
    Q_lsub_1[0] = Q_flux_LSUB
    # calculate p for LCC, add rcr term to pressure at LCC
    p_value_LCC = Q_flux_LCC*Rp_LCC
    temp1_lcc = 1/2*dt*(exp(-dt/(Rd_LCC*Cap_LCC))/Cap_LCC*Q_lcc_1[0] + Q_flux_LCC/Cap_LCC)
    temp2_lcc = exp(-dt/(Rd_LCC*Cap_LCC))*rcr_lcc[0]
    temp_lcc = temp1_lcc + temp2_lcc
    p_value_LCC += temp_lcc
    rcr_lcc[0] = temp_lcc
    Q_lcc_1[0] = Q_flux_LCC

    bcp_outflow_BCA = DirichletBC(Q, Constant(p_value_BCA), capid_markers, capid_seq['BCA'])
    bcp_outflow_LSUB = DirichletBC(Q, Constant(p_value_LSUB), capid_markers, capid_seq['LSUB'])
    bcp_outflow_LCC = DirichletBC(Q, Constant(p_value_LCC), capid_markers, capid_seq['LCC'])
    bcs['p']=[bcp_outflow, bcp_outflow_BCA, bcp_outflow_LCC, bcp_outflow_LSUB]

    #update variables
    NS_parameters['Q_bca_1'] = Q_bca_1
    NS_parameters['Q_lcc_1'] = Q_lcc_1
    NS_parameters['Q_lsub_1'] = Q_lsub_1
    NS_parameters['Q_aao_1'] = Q_aao_1
    NS_parameters['rcr_aao'] = rcr_aao
    NS_parameters['rcr_bca'] = rcr_bca
    NS_parameters['rcr_lcc'] = rcr_lcc
    NS_parameters['rcr_lsub'] = rcr_lsub
    globals().update(NS_parameters)

    #compute WSS
    if tstep % save_step == 0:
        sigma_ = sigma(u_,p_)
        T = -sigma_*n_normal
        Tn = inner(T,n_normal)
        Tt = T - Tn*n_normal
        Lt = inner(ww, Tt)*ds
        a_temp_wss = inner(shear_stress, ww) * ds
        A_wss = assemble(a_temp_wss, keep_diagonal=True)
        A_wss.ident_zeros()
        b_wss = assemble(Lt)
        solve(A_wss, shear_stress_.vector(), b_wss,'bicgstab','jacobi')
        # shear_stress_b_.interpolate(shear_stress_)

        

#define an intermediate function used for stabilization term
def dotun_(u,n):
    return (dot(u,n) - abs(dot(u,n)))/2
