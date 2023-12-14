#This script will take a surface mesh from a volumetric
#mesh and perform manipulation functions (e.g., coordinate
#changes, etc).

from utilities import *
import argparse
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class OasisSurfaceMeshManipulation():
	def __init__(self):
		#Read the volumetric mesh
		VolumeMesh=ReadVTUFile(MeshFileName)
		
		#Extract Surface Mesh
		MeshWall=ExtractSurface(ThresholdInBetweenCells(VolumeMesh,"CellEntityIds",1,1))
		
		#Compute Cell Normals
		MeshWall=SurfaceNormals(MeshWall)

	
	def ShrinkSurfacePoints(self,Surface,Tag="Shrink"):
		#Compute Mesh Quality to Get Size
		Surface=ComputeMeshQuality(Surface,"Size") #Return Cell Size 
		
		#Convert Cell data to Point Data
		CellToPoint = vtk.vtkCellDataToPointData()
		CellToPoint.SetInputData(Surface)
		CellToPoint.Update()
		Surface = CellToPoint.GetOutput()

		#Loop over all of the points
		for i in range(0,Surface.GetNumberOfPoints()):
			Factor_=Surface.GetPointData().GetArray("Quality").GetValue(i)
			if Tag=="Shrink": Factor_=Factor_*-1
			coord_=Surface.GetPoint(i)
			Nx_ =Surface.GetPointData().GetArray("Normals").GetValue(i*3)
			Ny_ =Surface.GetPointData().GetArray("Normals").GetValue(i*3+1)
			Nz_ =Surface.GetPointData().GetArray("Normals").GetValue(i*3+2)
			coord_new_=[coord_[0]+Nx_*Factor_,coord_[1]+Ny_*Factor_,coord_[2]+Nz_*Factor_]
			Surface.GetPoints().SetPoint(i,np.array(coord_new_))



if __name__=="__main__":
	OasisSurfaceMeshManipulation().Main() 
