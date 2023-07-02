#This code is written by Owais Khan on June 20th, 2023
#This code takes Oasis solution and inject contrast bolus

#from compute_flow_and_simulation_metrics import get_dataset_names
import json
from time import time
import argparse
from dolfin import *

class OasisAdvectionDiffusion():
	def __init__(self,Args):
		#Get all of the argumenets
		self.Args=Args
		#Number of Time Steps saved from the CFD Simulation
		self.Args.NumberOfTimesteps=int(((self.Args.EndStep-self.Args.StartStep)/self.Args.Increment))
		#Number of Time Steps to run Advection Diffusion Equation (5*CFD)
		self.Args.TotalTimesteps=int(self.Args.NumberOfTimesteps*self.Args.ContrastPeriodFactor)
		#Temporal Resolution of the CFD simulation
		self.Args.dt=self.Args.Period/self.Args.NumberOfTimesteps
		#Period for the Contrast Profile
		self.Args.PeriodContrast=self.Args.Period*self.Args.ContrastPeriodFactor

		#print the results
		if MPI.rank(MPI.comm_world) == 0:
			print ("-"*75)
			print ("\n")
			print ("Mesh File:              %s"%self.Args.MeshFileName)
			print ("Velocity File           %s"%self.Args.InputFileName)
			print ("\n")
			print ("Diffusion Coeff:        %.05f"%self.Args.DiffusionCoefficient)
			print ("Period:                 %.05f"%self.Args.Period)
			print ("Contrast Period Factor: %.05f"%self.Args.ContrastPeriodFactor)
			print ("No of Timesteps:        %.05f"%self.Args.NumberOfTimesteps)
			print ("Velocity Order:         %d"%self.Args.VelocityOrder)
			print ("\n")
			print ("Start TimeStep for CFD:  %s"%self.Args.StartStep)
			print ("Ending TimeStep for CFD: %s"%self.Args.EndStep)
			print ("Increment Step:          %s"%self.Args.Increment)
			
	def Main(self):
		#--------------------------------------------------------
		print ("-"*75)
		#Define all of the input parameters
		#Mesh
		print ("--- Loading Mesh File: %s"%self.Args.MeshFileName)
		self.Mesh=Mesh(self.Args.MeshFileName) #Mesh
		
		#Temporal Resolutions
		print ("--- Creating Temporal Resolution Variable.")
		self.Args.dt=self.Args.Period/self.Args.NumberOfTimesteps #Temp Res
		k=Constant(self.Args.dt)

		#Assign a diffusion constant
		print ("--- Creating Diffusivity Coefficient Variable.")
		D=Constant(self.Args.DiffusionCoefficient)		
		
		#--------------------------------------------------------
		print ("-"*75)
		#Create Function Space and Test/Trial Functions
		print ("--- Creating a Vector Function Scape for Velocity.")
		W=VectorFunctionSpace(self.Mesh,"CG",self.Args.VelocityOrder) # For velocity
		print ("--- Creating a Scalar Function Space for Contrast Concentration.")
		V=FunctionSpace(self.Mesh,"CG",2) #For Contrast


		print ("--- Creating FEM Test Function.")
		v= TestFunction(V) #Test Function
	
		print ("--- Creating Velocity Function.")	
		w  = Function(W) #Velocity Function
		print ("--- Creating Contrast Concentration Function.")
		u  = TrialFunction(V) #Concentration Function
		#u_n=Function(V)
		
		#---------------------------------------------------------
		print ("-"*75)
		#Assing Boundary Conditions
		print ("--- Creating Boundary Mesh Function.") 
		#Read the Boundary Mesh
		Boundary=MeshFunction("size_t", self.Mesh, self.Mesh.geometry().dim()-1, self.Mesh.domains())

		print ("--- Reading the Mesh Info file for Inlet/Outlet Ids.") 
                # Read Case Parameters from Simulation File
		info = self.Args.MeshFileName.split(".xml")[0] + "_info.json"
		with open(info) as f: info = json.load(f)
		id_in=[]
		id_in[:] = info['inlet_id']
		id_out=[]
		id_out[:] = info['outlet_ids']
 		
		print ("-"*75)
		#Define time-dependent Contrast Boundary Condition #Eslami+, JBioMechEng, 2019
		print ("--- Creating Expression for Contrast Concentration Profile.")
		ConcentrationEquation = Expression("cmin+0.5*(cmax-cmin)*(1-cos(pi*((t-Ts)/Td)))",cmin=0.0,cmax=1.0,t=0.0,Ts=0.0,Td=self.Args.PeriodContrast,degree=2)
		
		#Assing the inflow boundary for the concentration
		print ("--- Assigning Wall Boundary Condition.")
		Wall = DirichletBC(V, Constant(0.0), Boundary,0) #Assuming wall=0
		print ("--- Assigning Dirichlet Boundary Condition. Assume InletId=1")
		Inflow = DirichletBC(V,ConcentrationEquation,Boundary,1) #Assuming inflow=1
		BoundaryConditions = [Wall,Inflow]
		u_n = interpolate(ConcentrationEquation,V)
		print ("\n")

		#---------------------------------------------------------
		#Write the Function to minimize
		print ("--- Creating Variational Equation")
		F = ((u - u_n) / k)*v*dx + dot(w, grad(u))*v*dx	- D*dot(grad(u), grad(v))*dx 
		print ("--- Separating LHS and RHS")
		a1=lhs(F)
		L1=rhs(F)
		A1=assemble(a1)
		u = Function(V)
		#Load the Velocity File Series from CFD simulation
		f_u = HDF5File(MPI.comm_world, self.Args.InputFileName, "r")
		#f_u.read(w,f"/velocity/vector_{self.Args.NumberOfTimesteps}")
		#----------------------------------------------------------
		print ("-"*75)
		# Time-stepping
		t = 0.0
		for i in range(self.Args.TotalTimesteps):
			# Update current time
			Progress_=(t/(self.Args.Period*self.Args.ContrastPeriodFactor))*100
			print ("Current Time is: %.05f. Completed: %s"%(t,Progress_))
			t += self.Args.dt
   
			#Update the Contrast Inflow 
			ConcentrationEquation.t=t	
 
			# Read velocity from file
			f_u.read(w,f"/velocity/vector_{self.Args.NumberOfTimesteps}")

			# Solve variational problem for time step
			b1=assemble(L1)
			[bc.apply(A1,b1) for bc in BoundaryConditions]
			#solve(a1 == L1, u, BoundaryConditions)#A1,c.vector(),b1,"gmres","default")
			solve(A1,c.vector(),b1)#,"gmres","default")
			# Update previous solution
			u_n.assign(u)


if __name__=="__main__":
	parser = argparse.ArgumentParser(description="This script will run an advection-diffusion equation on Oasis generated velocity field.")

	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The file containing the velocity data generated by oasis.")
	
	parser.add_argument('-Period', '--Period', type=float, required=False, default=1.0, dest="Period",help="The duration of one cardiac cycle.")
	
	parser.add_argument('-ContrastPeriodFactor', '--ContrastPeriodFactor', type=float, required=False, default=5.0, dest="ContrastPeriodFactor",help="The duration of the contrast bolus injection to peak contrast as a multiple of Period. Default is 5 (i.e., 5*Period)")
        
	parser.add_argument('-DiffusionCoefficient', '--DiffusionCoefficient', type=float, required=False, default=0.04, dest="DiffusionCoefficient",help="The diffusivity coefficient. Default is 0.04, which assumes Schmit Number (nu/D) of 1.")

	parser.add_argument('-VelocityOrder', '--VelocityOrder', type=int, required=False, default=1, dest="VelocityOrder",help="The polynomial order for the velocity field.")
	
	parser.add_argument('-StartStep', '--StartStep', type=int, required=False, default=0, dest="StartStep",help="The starting timestep from the CFD simulation. Default=20000")
	
	parser.add_argument('-EndStep', '--EndStep', type=int, required=False, default=1000, dest="EndStep",help="The ending timestep for CFD simulation.")
	
	parser.add_argument('-Increment', '--Increment', type=int, required=False, default=1, dest="Increment",help="The Increment Between Timestep Saving.")
	
	parser.add_argument('-MeshFileName', '--MeshFileName', type=str, required=True, dest="MeshFileName",help="The filename that contains the mesh file")
        
	args=parser.parse_args()
	OasisAdvectionDiffusion(args).Main()
	
