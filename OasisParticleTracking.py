from dolfin import *
from utilities import *
import argparse
import vtk

class OasisParticleTracking():
	def __init__(self,args):
		self.Args=args
		self.Args.dt=self.Args.Period/(self.Args.Stop-self.Args.Start)
	
	def Main(self):
		#Read the mesh
		Mesh=ReadFenicsMesh(self.Args.InputMeshFile)
		print ("Mesh Loaded...")
	
		#Read the velocity data	
		f = HDF5File(MPI.comm_world, self.Args.InputVelocityFile, "r")
                # Get names of data to extract
		Start = 0
		if MPI.rank(MPI.comm_world) == 0:
			print("The post processing Starts from", Start)
		dataset_names = get_dataset_names(f, Start=self.Args.Start,num_files=self.Args.Stop)
    
		# Create Function space
		print ("Creating Function Spaces")
    		V = VectorFunctionSpace(mesh, "CG", velocity_degree)
		u = Function(V)
	
		#Loop over the velocity field
		dt=self.Args.dt
		for dataset_name_ in dataset_names:
			f.read(u,dataset_name_)
			



if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will seed points near the given surface and track them until they reach the outlets.")

	parser.add_argument('-InputVelocityFile', '--InputVelocityFile', type=str, required=True, dest="InputVelocityFile",help="The h5 file containing the velocity field.")
        
	parser.add_argument('-VelocityOrder', '--VelocityOrder', type=int, required=False, default=1, dest="VelocityOrder",help="The polynomial order for the velocity field.")
        
	parser.add_argument('-InputMeshFile', '--InputMeshFile', type=str, required=True, dest="InputMeshFile",help="The file contains the mesh in vtu format.")
	
	parser.add_argument('-Start', '--Start', type=int, required=False, dest="Start",default=1,help="The Starting time step")
	
	parser.add_argument('-Stop', '--Stop', type=int, required=False, dest="Stop",default=1000,help="The ending time step")
	
	parser.add_argument('-Period', '--Period', type=float, required=False, dest="Perioid",default=0.951,help="The length of the cardiac cycle.")
	
 
        #Output Filename 
	parser.add_argument('-OutputFile', '--OutputFile', type=str, required=False, dest="OutputFile",help="The output file in which to store the dolfin mesh.")

                
	args=parser.parse_args()
        
	OasisParticleTracking(args).Main()


