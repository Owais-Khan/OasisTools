import os, sys
import vtk, numpy, h5py
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
import argparse
from glob import glob

class Hdf5ToVtu():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.VelocityFile is None:
			self.Args.VelocityFile=glob(self.Args.InputFolder+"/*.h5")[0]
		if self.Args.MeshFile is None:
			self.Args.MeshFile=self.Args.InputFolder+"/mesh.vtu"
	def Main(self):
		
		#Read the VTU File
		print ("--- Reading %s:"%self.Args.MeshFile)
		Mesh=self.ReadVTUFile(self.Args.MeshFile)
		
		#Read the H5 file 
		print ("--- Reading %s:"%self.Args.VelocityFile)
		Velocity,Pressure=self.ReadH5(self.Args.VelocityFile)
	
		#Project Velocity to VTU File
		print ("--- Convert to VTU")	
		Mesh=self.ProjectVelocityToMesh(Mesh,Velocity,"Velocity")
		Mesh=self.ProjectVelocityToMesh(Mesh,Pressure,"Pressure")

		#Save the File
		print ("--- Saving VTU File: %s"%self.Args.VelocityFile.replace(".h5",".vtu"))
		self.WriteVTUFile(self.Args.VelocityFile.replace(".h5",".vtu"),Mesh)

	def ProjectVelocityToMesh(self,Mesh,Data,Name):
		DataArray=numpy_to_vtk(Data)
		DataArray.SetName(Name)
		Mesh.GetPointData().AddArray(DataArray)
		return Mesh

	def ReadH5(self,FileName):
		f=h5py.File(FileName,'r')
		a_group_key = list(f.keys())[0]
		keys = list(f[a_group_key])
		velocity=f[a_group_key]['u'][:]
		pressure=f[a_group_key]['p'][:]
		return velocity,pressure	

	def WriteVTUFile(self,FileName,Data):
		writer=vtk.vtkXMLUnstructuredGridWriter()
		writer.SetFileName(FileName)
		writer.SetInputData(Data)
		writer.Update()

	def ReadVTUFile(self,FileName):
		reader=vtk.vtkXMLUnstructuredGridReader()
		reader.SetFileName(FileName)
		reader.Update()
		return reader.GetOutput()

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="Convert HD5F written by Mehdi back to Vtu Format and extract random probe points.")

	parser.add_argument('-MeshFile', '--MeshFile', type=str,required=False, dest="MeshFile",help="Filename that contains the mesh in vtu format ")

	parser.add_argument('-VelocityFile', '--VelocityFileName', type=str, required=False, dest="VelocityFile", help="The file that contains velocity.")
	
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=False, dest="InputFolder", help="The input folder contains h5 and vtu.")

	args=parser.parse_args()
	Hdf5ToVtu(args).Main()
