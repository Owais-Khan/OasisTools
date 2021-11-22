import dolfin as df
import numpy as np
import argparse
from glob import glob
import h5py

mpi_comm = df.MPI.comm_world
my_rank = df.MPI.rank(mpi_comm)


class Utilities():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Convert h5 to vtu
		if self.Args.ConvertHDF5ToVTU==True:
			#Loop over all files if no start/stop/increment provided
			if self.Args.StartTimestep is None:
				FileNames_=sorted(glob(self.Args.InputFolder+"/*up.h5"))
				self.Args.StartTimestep=int(FileNames_[0].split("=")[1].split("_")[0])
				self.Args.StopTimestep=int(FileNames_[-1].split("=")[1].split("_")[0])
				self.Args.Increment=int(FileNames_[1].split("=")[1].split("_")[0])-self.Args.StartTimestep
				print ("--- The Start Timestep is: %d"%self.Args.StartTimestep)
				print ("--- The End Timestep is:   %d"%self.Args.StopTimestep)
				print ("--- The Increment is:      %d"%self.Args.Increment)
				
			print ("--- Converting Velocity from HDF5 to VTU format")
			self.ConvertHDF5ToVTU()
		
		if self.Args.ProjectVelocityToNewMesh==True:
			print ("--- Projecting Velocity on New Mesh")
			self.ProjectVelocityToNewMesh()
			
	
	#This function will convert HDF5 file format to VTU File Format
	def ConvertHDF5ToVTU(self):
		#Load the Mesh
		print ("------ Loading Mesh %s"%self.Args.MeshFileName)
		Mesh=df.Mesh(self.Args.MeshFileName)
		
		#Create a Function Space
		V=df.VectorFunctionSpace(Mesh,"CG",self.Args.VelocityOrder)
		u=df.Function(V)
		Npts=u.vector().size()/3

		#Loop over the time steps
		for i in range(self.Args.StartTimestep,self.Args.StopTimestep+self.Args.Increment,self.Args.Increment):
			FileName_=glob(self.Args.InputFolder+"%.05d*.h5")[0]
			print ("------ Looping over: %s"%FileName_)
			Velocity_=df.HDF5File(Mesh.MPI.comm_world,FileName_, "r")
			Velocity_.read(u,"u")		


	def ProjectVelocityToNewMesh(self):
		#Load the Reference Mesh
		print ("------ Loading Mesh: %s"%self.Args.MeshFileName)
		Mesh=df.Mesh(self.Args.MeshFileName)

		#Load New Mesh
		VV=VectorFunctionSpace(mesh,"CG",self.uOrder)
		u=Function(VV)
		N=u.vector().size()/3	
		






if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This is a general Utilitiy script to perform several taks related to Oasis solver, pre-processing and post-processing.")
	
	parser.add_argument('-MeshFileName', '--MeshFileName', type=str,required=True, dest="MeshFileName",help="The path/to/filename where ")
	
	parser.add_argument('-VelocityOrder', '--VelocityOrder', type=int,required=True, dest="VelocityOrder",help="The polynomial order for the velocity function space.")
	
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The folder that containts the velocity data")


#---------------------------- Convert HDF5 Files to VTU Files ------------        
	parser.add_argument('-ConvertHDF5ToVTU', '--ConvertHDF5ToVTU', type=bool, required=False, default=False, dest="ConvertHDF5ToVTU",help="Convert velocity from H5 format to Vtu format")
	
	parser.add_argument('-StartTimestep', '--StartTimestep', type=int, required=False, dest="StartTimestep",help="Starting timestep to process the data")
	parser.add_argument('-StopTimestep', '--StopTimestep', type=int, required=False, dest="StopTimestep",help="Last timestep to process the data")
	parser.add_argument('-Increment', '--Increment', type=int, required=False, dest="Increment",help="Timestep increment")





#----------------------------- ProjectVelocityToNewMesh------------------


	#parser.add_argument('-ProjectVelocityToNewMesh', '--ProjectVelocityToNewMesh', type=bool, required=False, default=True, dest="ProjectVelocityToNewMesh",help="Project the Velocity and Pressure from one mesh to another mesh.")

	#parser.add_argument('-MeshFileName', '--MeshFileName', type=str,required=False, dest="MeshFileName",help="The path/to/filename where ")

	#parser.add_argument('-NewMeshFileName', '--NewMeshFileName', type=str, required=False, dest="NewMeshFileName", help="Related to ProjectVelocityToNewMesh. Filename for the new mesh on which to project the velocity data")
	
	#parser.add_argument('-OutFolder', '--OutFolder', type=str, required=False, dest="OutFolder", help="The output folder to store all the projected velocity data")

#------------------------------- Velocity Data -------------------------
	#parser.add_argument('-VelocityFolder', '--VelocityFolder', type=str, required=False, dest="VelocityFolder", help="The folder where velocity data is stored")

	#parser.add_argument('-StartTimeStep', '--StartTimeStep', type=int, required=False, dest="StartTimeStep", help="The starting time step to be processed")
	
	#parser.add_argument('-EndTimeStep', '--EndTimeStep', type=int, required=False, dest="EndTimeStep", help="The end time step to be processed")
	
	#parser.add_argument('-Increment', '--Increment', type=int, required=False, dest="Increment", help="Timestep increment")

	args=parser.parse_args()
	Utilities(args).Main()
