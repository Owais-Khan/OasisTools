import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from utilities import *

class OasisMeshWriterForSimVascular():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Lead the Vtu file from mesh-complete folder
		MeshVTK=ReadVTUFile(self.Args.InputFolder+"/mesh-complete.mesh.vtu")
		WallVTK=ReadVTPFile(self.Args.InputFolder+"/walls_combined.vtp")		

		#Read all of the inlet and outlet ids and remove wall files
		BoundaryFilesAll=sorted(glob(self.Args.InputFolder+"/mesh-surfaces/*.vtp"))
		BoundaryFiles=[]
		
		#Append the inflow.vtp to the start of the list
		try:
			idx=BoundaryFilesAll.index(self.Args.InputFolder+"/mesh-surfaces/inflow.vtp")
			BoundaryFiles.append(BoundaryFilesAll[idx])
			BoundaryFilesAll.remove(BoundaryFiles[0])
		except:
			print ("--- No inflow.vtp file found. Can't decide on inflow")
			print ("--- Exiting...")
			exit(1)
			
		#Append all inlets/outlets except the wall	
		for BoundaryFile_ in BoundaryFilesAll:
			if BoundaryFile_.find("wall")<0: BoundaryFiles.append(BoundaryFile_)

		#Store all of the CellIds of Inlet and Outlet
		Array=np.zeros(MeshVTK.GetNumberOfCells()) #0 ====> Tetrahedral Cells inside the volume
		
		print ("--- Assigning Wall Element ID")
		for i in range(0,WallVTK.GetNumberOfCells()):
			CellID_ =WallVTK.GetCellData().GetArray("GlobalElementID").GetValue(j)-1
			





	
		print ("--- Assigning Inlet and Outlet Files in mesh-complete/mesh-surfaces/")
		for i in range(2,len(BoundaryFiles)+2):
			BoundarySurface_=ReadVTPFile(BoundaryFiles[i-2])
			NCells_=BoundarySurface_.GetNumberOfCells()
			print ("------ No. of Cells in %s are: %d"%(BoundaryFiles[i-2],NCells_))
			for j in range(NCells_):
				ID_=BoundarySurface_.GetCellData().GetArray("GlobalElementID").GetValue(j)-1
				Array[ID_]=i


		#Add the CellEntityIds array to MeshVTK
		CellEntityIds=numpy_to_vtk(Array)
		CellEntityIds.SetName("CellEntityIds")
		MeshVTK.GetCellData().AddArray(CellEntityIds)	

		#Add a new array for CellEntityIds
		MeshVTK.GetCellData().RemoveArray("GlobalElementID")		
		MeshVTK.GetCellData().RemoveArray("ModelRegionID")		
		MeshVTK.GetPointData().RemoveArray("GlobalNodeID")
		
		#Write Output File

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take a mesh-complete folder from SimVascular and write a dolfin mesh file. The CellEntityIds 0 is for volume, 1 is for wall, 2 is for inlet and remaining is for outlets.")
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The path to mesh-complete folder from SimVascular")
        
	#Output Filename 
	parser.add_argument('-OutputFile', '--OutputFile', type=str, required=False, dest="OutputFile",help="The output file in which to store the dolfin mesh.")

	args=parser.parse_args()
	OasisMeshWriterForSimVascular(args).Main()



 
