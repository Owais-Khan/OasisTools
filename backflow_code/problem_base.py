__author__ = "Yiyang Fu <yiyang.fu@mail.utoronto.ca>"
__date__ = "2022-04-14"

# from dolfin import Mesh, MeshFunction
# import csv
# from scipy.optimize import curve_fit
# from math import acos, sqrt
from dolfin import *
from numpy import pi
from scipy.interpolate import interp1d
import numpy as np
# import mpi4py


# Read mesh from xml file
mesh = Mesh()
mesh_filename = 'case11_ab_final_u2_mesh.xdmf'
with XDMFFile(mesh_filename) as XDMFoutfile:
	XDMFoutfile.read(mesh)
	XDMFoutfile.close()

boundaries_filename = 'case11_ab_final_u2_mesh_facet_markers.xdmf'
capid_markers_model = MeshValueCollection('size_t',mesh,mesh.topology().dim()-1)
with XDMFFile(boundaries_filename) as XDMFoutfile:
	XDMFoutfile.read(capid_markers_model)
	XDMFoutfile.close()

capid_markers = cpp.mesh.MeshFunctionSizet(mesh, capid_markers_model)

capid_seq = {
	'inlet': 6,
	'outlet': 7,
	'BCA': 9,
	'LCC': 8,
	'LSUB': 5
}


ds = Measure('ds', domain=mesh, subdomain_data=capid_markers)
# inlet_area=assemble(Constant(1)*ds(capid_seq['inlet']))
outlet_area_Aao=assemble(Constant(1)*ds(capid_seq['outlet']))
outlet_area_BCA=assemble(Constant(1)*ds(capid_seq['BCA']))
outlet_area_LSUB=assemble(Constant(1)*ds(capid_seq['LSUB']))
outlet_area_LCC=assemble(Constant(1)*ds(capid_seq['LCC']))
#RCR Model Parameters
R_tot = 1500.0
Area_total = outlet_area_Aao + outlet_area_BCA + outlet_area_LSUB + outlet_area_LCC
R_value = R_tot*Area_total/outlet_area_Aao
R_value_BCA = R_tot*Area_total/outlet_area_BCA
R_value_LSUB = R_tot*Area_total/outlet_area_LSUB
R_value_LCC = R_tot*Area_total/outlet_area_LCC


cap_tot = 0.001
Rp_BCA = R_value_BCA*0.09
Cap_BCA = cap_tot*outlet_area_BCA/Area_total
Rd_BCA = R_value_BCA*0.91

Rp_LCC = R_value_LCC*0.09
Cap_LCC = cap_tot*outlet_area_LCC/Area_total
Rd_LCC = R_value_LCC*0.91

Rp_LSUB = R_value_LSUB*0.09
Cap_LSUB = cap_tot*outlet_area_LSUB/Area_total
Rd_LSUB = R_value_LSUB*0.91

Rp_aao = R_value*0.09
Cap_aao = cap_tot*outlet_area_Aao/Area_total
Rd_aao = R_value*0.91

Q_bca_1 = [0.0]
Q_lcc_1 = [0.0]
Q_lsub_1 = [0.0]
Q_aao_1 = [0.0]

rcr_aao = [0.0]
rcr_bca = [0.0]
rcr_lsub = [0.0]
rcr_lcc = [0.0]


mu = 0.037736 # gr*cm^-1*s^-1
vector = VectorFunctionSpace(mesh,"CG", 1)
ww = TestFunction(vector)
shear_stress = TrialFunction(vector)

shear_stress_ = Function(vector)
shear_stress_.rename('shear_stress','shear_stress')
shear_stress_.set_allow_extrapolation(True)

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))



