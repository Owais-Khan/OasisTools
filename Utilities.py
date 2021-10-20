import dolfin as df
import numpy as np
import argparse

mpi_comm = df.MPI.comm_world
my_rank = df.MPI.rank(mpi_comm)


class Utilities():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		if self.Args.ProjectVelocityToNewMesh==True:
			print ("--- Projecting Velocity on New Mesh")
			self.ProjectVelocityToNewMesh()

	def ProjectVelocityToNewMesh(self):
		#Load the Reference Mesh
		print ("------ Loading Mesh: %s"%self.Args.MeshFileName)
		Mesh=df.Mesh(self.Args.MeshFileName)

		#Load New Mesh

	def 



	def LoadMesh(self,FileName):
		print ("Loading Mesh: %s"%FileName)
		Mesh = df.Mesh(FileName)
		return Mesh	
	







if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This is a general Utilitiy script to perform several taks related to Oasis solver, pre-processing and post-processing.")

#----------------------------- ProjectVelocityToNewMesh------------------


	parser.add_argument('-ProjectVelocityToNewMesh', '--ProjectVelocityToNewMesh', type=bool, required=False, default=True, dest="ProjectVelocityToNewMesh",help="Project the Velocity and Pressure from one mesh to another mesh.")

	parser.add_argument('-MeshFileName', '--MeshFileName', type=str,required=False, dest="MeshFileName",help="The path/to/filename where ")

	parser.add_argument('-NewMeshFileName', '--NewMeshFileName', type=str, required=False, dest="NewMeshFileName", help="Related to ProjectVelocityToNewMesh. Filename for the new mesh on which to project the velocity data")
	
	parser.add_argument('-OutFolder', '--OutFolder', type=str, required=False, dest="OutFolder", help="The output folder to store all the projected velocity data")

#------------------------------- Velocity Data -------------------------
	#parser.add_argument('-VelocityFolder', '--VelocityFolder', type=str, required=False, dest="VelocityFolder", help="The folder where velocity data is stored")

	#parser.add_argument('-StartTimeStep', '--StartTimeStep', type=int, required=False, dest="StartTimeStep", help="The starting time step to be processed")
	
	#parser.add_argument('-EndTimeStep', '--EndTimeStep', type=int, required=False, dest="EndTimeStep", help="The end time step to be processed")
	
	#parser.add_argument('-Increment', '--Increment', type=int, required=False, dest="Increment", help="Timestep increment")

	args=parser.parse_args()
	Utilities(args).Main()
