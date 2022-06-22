import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from utilities import *
from dolfin import *
from ufl.tensors import ListTensor


parameters["reorder_dofs_serial"] = False

class OasisConvertSimVascularResults():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			os.system("mkdir %s/../Results_OasisFormat"%self.Args.InputFolder)
			self.Args.OutputFolder="%s/../Results_OasisFormat"%self.Args.InputFolder

	def Main(self):
		#Read all of the filenames
		FileNames=sorted(glob(self.Args.InputFolder+"/*.vtu"))
		print ("--- Found %d files in %s"%(len(FileNames),self.Args.InputFolder))

		#Load the dolfin mesh
		print ("--- Mesh loading from %s"%self.Args.MeshFileName)
		Mesh_=Mesh(self.Args.MeshFileName)

		print ("--- Writing mesh to: %s"%self.Args.OutputFolder)
		File("%s/mesh.xml.gz"%self.Args.OutputFolder)<<Mesh_

		#Create a function space
		print ("--- Creating function spaces")
		VV=VectorFunctionSpace(Mesh_,"CG",1)
		V=FunctionSpace(Mesh_,"CG",1)
		u=Function(VV)
		q=Function(V)

		#Create an HDF5 File to store the velocity
		VelocityFile_h5 = HDF5File(MPI.comm_world,self.Args.OutputFolder+"/u.h5", 'w')
		PressureFile_h5 = HDF5File(MPI.comm_world,self.Args.OutputFolder+"/p.h5", 'w')
		VelocityMean_h5 = HDF5File(MPI.comm_world,self.Args.OutputFolder+"/u_mean.h5", 'w')

		#Compute an array for average velocity
		u_mean=0
		
		#Loop over all of the files
		for i in range(len(FileNames)):
			print ("--- Looping over: %s"%FileNames[i])

			#Load the pressure and velocity data
			Data_=ReadVTUFile(FileNames[i])
			velocity_=vtk_to_numpy(Data_.GetPointData().GetArray("velocity"))
			pressure_=vtk_to_numpy(Data_.GetPointData().GetArray("pressure"))

			#Compute the mean velocity
			u.vector()[:]=velocity_.flatten("F")
			q.vector()[:]=pressure_.flatten("F")

			if i==0:
				u_mean=velocity_.flatten("F")
			else:
				u_mean+=velocity_.flatten("F")
			
			#Append it to the velocity and pressure file
			VelocityFile_h5.write(u,"velocity",float(i))
			PressureFile_h5.write(q,"pressure",float(i))
		
		VelocityMean_h5.write(u,"u_mean",float(i))
		
		del VelocityFile_h5
		del PressureFile_h5
		del VelocityMean_h5
		
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take results from SimVascular and convert into Oasis format.")
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The path to folder containing the results in .vtu format")
	parser.add_argument('-MeshFileName', '--MeshFileName', type=str, required=True, dest="MeshFileName",help="The file containing the mesh in .xml.gz format. Use vmtkmeshwriter to convert any results file into a mesh file.")
       
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder to store the converted velocity and pressures (i.e., u.h5 and p.h5).")

	args=parser.parse_args()
	OasisConvertSimVascularResults(args).Main()
