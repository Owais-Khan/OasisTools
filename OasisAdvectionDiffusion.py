import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from compute_flow_and_simulation_metrics import get_dataset_names
import numpy as np
import vtk
import argparse
from utilities import *
from dolfin import *
from ufl.tensors import ListTensor

class OasisAdvectionDiffusion():
	def __init__(self,Args):
		#Get all of the argumenets
		self.Args=Args
		
		#print the results
		if MPI.rank(MPI.comm_world) == 0:
			print ("Mesh File:        %s"%self.Args.MeshFileName)
			print ("Velocity File     %s"%self.Args.VelocityFileName) #Contains Velocity from Oasis
			print ("Diffusion Coeff:  %.05f"%self.Args.DiffusionCoefficient)
			print ("Period:           %.05f"%self.Args.Period)
			print ("No of Timesteps:  %.05f"%self.Args.NumerOfTimesteps)
			print ("Velocity Order:   %d"%self.Args.VelocityOrder)
			print ("Number of Cycles: %d"%self.Args.NumberOfCycles) 
	def Main(self):
		#Load the Mehs
		mesh=Mesh(self.Args.InputFileName)

		#Define the temporal resolution
		self.Args.dt=self.Args.Period/self.Args.NumberOfTimesteps

		#Get the temporal resolution
		k=Constant(self.Args.dt)

		#Assign a diffusion constant
		D=Constant(Self.Args.DiffusionCoefficient)		
		
		#Create Function Space and Test/Trial Functions
		VV=VectorFunctionSpace(mesh,"CG",self.Args.VelocityOrder) #Velocity Function Space
		VC=FunctionSpace(mesh,"CG",2) #Contrast Function Space
		
		u  =Function(VV) #Velocity Function
		c  =Function(VC) #Concentration Function
		c_n=Function(VC)
		v=Function(VC) #Test Function

		#Write the Function to minimize
		F = ((c - c_n) / k)*v*dx + dot(u, grad(c))*v*dx	+ D*dot(grad(c), grad(v))*dx 

		#Load the Velocity File Series
		f = HDF5File(MPI.comm_world, self.Args.VelocityFileName, "r")
		dataset_names = get_dataset_names(f, start=start)
		
		# Time-stepping
		t = 0
		for n in range(num_steps):
			# Update current time
			t += dt
    
			# Read velocity from file
			timeseries_w.retrieve(w.vector(), t)

			# Solve variational problem for time step
			solve(F == 0, u)

			# Save solution to file (VTK)

    
			# Update previous solution
			u_n.assign(u)

    
			# Update progress bar
			progress.update(t / T)


			
