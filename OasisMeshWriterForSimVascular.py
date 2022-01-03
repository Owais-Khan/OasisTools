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
		#SurfaceVTK=ReadVTPFile(self.Args.InputFolder+"/mesh-complete.exterior.vtp")		

		#Read all of the inlet and outlet ids and remove wall files
		BoundaryFilesAll=sorted(glob(self.Args.InputFolder+"/mesh-surfaces/*.vtp"))
		BoundaryFiles=[]
		try: #Append the inflow.vtp to the start of the list
			idx=BoundaryFilesAll.index(self.Args.InputFolder+"/mesh-surfaces/inflow.vtp")
			BoundaryFiles.append(BoundaryFilesAll[idx])
			BoundaryFilesAll.remove(BoundaryFiles[0])
		except:
			print ("--- No inflow.vtp file found. Can't decide on inflow")
			print ("--- Exiting...")
			exit(1)
                
		#Read All of the inlet and outlet files
		BoundaryFiles=[self.Args.InputFolder+"/walls_combined.vtp",self.Args.InputFolder+"/mesh-surfaces/inflow.vtp"] 
		for BoundaryFile_ in BoundaryFilesAll:
			if BoundaryFile_.find("wall")<0:
				BoundaryFiles.append(BoundaryFile_)

		#Loop over all of the surface files and read the surface node ids
		BoundaryNodeIds=[]
		BoundaryEntityIds=[]
		counter=1
		for BoundaryFile_ in BoundaryFiles:
			BoundarySurface_=ReadVTPFile(BoundaryFile_)
			BoundaryNodeIds+=[BoundarySurface_.GetPointData().GetArray("GlobalNodeID").GetValue(i) for i in range(BoundarySurface_.GetNumberOfPoints())]
			BoundaryEntityIds+=[counter]*BoundarySurface_.GetNumberOfPoints()
			counter+=1 

		BoundaryNodeIds=np.array(BoundaryNodeIds)
		BoundaryEntityIds=np.array(BoundaryEntityIds)

		#Create Cell Entity IDs
		CellEntityIds=[0]*MeshVTK.GetNumberOfCells()
	
		#Create a triangle cell for the boundary
		triangle_=vtk.vtkTriangle()
		counter=1
		for BoundaryFile_ in BoundaryFiles:
			SurfaceVTK=ReadVTPFile(BoundaryFile_)
			for i in range(SurfaceVTK.GetNumberOfCells()):
				#Get the point Ids for the triangle
				pointIds_=[SurfaceVTK.GetCell(i).GetPointIds().GetId(j) for j in range(3)]
				pointIdsGlob_=[SurfaceVTK.GetPointData().GetArray("GlobalNodeID").GetValue(j)-1 for j in pointIds_]
				#Set the triangle ids
				triangle_.GetPointIds().SetId(0,pointIdsGlob_[0])
				triangle_.GetPointIds().SetId(1,pointIdsGlob_[1])
				triangle_.GetPointIds().SetId(2,pointIdsGlob_[2])
	
				#Add to the MeshVTK as Next Cell
				MeshVTK.GetCells().InsertNextCell(triangle_)

             			#Add the Mesh Type to the Array
				MeshVTK.GetCellTypesArray().InsertNextValue(5)
		
				#Add the values of the Array
				MeshVTK.GetCellData().GetArray("GlobalElementID").InsertNextValue(SurfaceVTK.GetCellData().GetArray("GlobalElementID").GetValue(i))
				
				#Add CellEntityIds
				CellEntityIds.append(counter)
			counter+=1

		#Add CellEntityIds
		CellEntityIds=np.array(CellEntityIds)
		CellEntityIdsVTK=numpy_to_vtk(CellEntityIds)
		CellEntityIdsVTK.SetName("CellEntityIds")
		MeshVTK.GetCellData().AddArray(CellEntityIdsVTK)
		
		#Loop over all of the Boundary Files 
		MeshVTK.GetCellData().RemoveArray("ModelRegionID")
		MeshVTK.GetCellData().RemoveArray("GlobalElementID")
		MeshVTK.GetPointData().RemoveArray("GlobalNodeID")

		#Write the Mesh File in VTK and XML format
		WriteVTUFile(self.Args.InputFolder+".vtu",MeshVTK)

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take a mesh-complete folder from SimVascular and write a dolfin mesh file. The CellEntityIds 0 is for volume, 1 is for wall, 2 is for inlet and remaining is for outlets.")
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The path to mesh-complete folder from SimVascular")
        
	#Output Filename 
	parser.add_argument('-OutputFile', '--OutputFile', type=str, required=False, dest="OutputFile",help="The output file in which to store the dolfin mesh.")

	args=parser.parse_args()
	OasisMeshWriterForSimVascular(args).Main()



 
