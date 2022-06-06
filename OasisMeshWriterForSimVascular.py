#This script was written by Owais Khan on January 3, 2022 
#Cardiovascular Imaging, Modeling and Biomechanics Lab (CIMBL)
#Ryerson University, Canada

import json
import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from Utilities import *
from vmtk import vtkvmtk, vmtkscripts

class OasisMeshWriterForSimVascular():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFileName is None:
			self.Args.OutputFileName=self.Args.InputFolder

		os.system("rm %s.xml.gz"%self.Args.OutputFileName)
		os.system("rm %s.vtu"%self.Args.OutputFileName)
		os.system("rm %s_info.json"%self.Args.OutputFileName)

		self.IDs={"check_surface": True}

	def Main(self):
		if self.Args.ScalingFactor is not None:
			#Get All of the files in the mesh complete folder
			print ("--- Scaling the mesh from cm to mm since SimVascular is usually in cm")
			print ("--- Using the Scaling Factor: %.02f"%self.Args.ScalingFactor)
			AllFilesVTP=glob("%s/*.vtp"%self.Args.InputFolder)+glob("%s/mesh-surfaces/*.vtp"%self.Args.InputFolder)
			for File_ in AllFilesVTP:
				print ("------ Scaling %s"%(File_))
				os.system("vmtksurfacescaling -ifile  %s -ofile %s -scale %f"%(File_,File_,self.Args.ScalingFactor))
			os.system("vmtkmeshscaling -ifile %s/mesh-complete.mesh.vtu -ofile %s/mesh-complete.mesh.vtu -scale %f"%(self.Args.InputFolder,self.Args.InputFolder,self.Args.ScalingFactor))
		else:
			print ("------------- NOT SCALING THE MESH-COMPLETE FOLDER --------------")
			print ("---------- MAKE SURE THE MESH-COMPLETE is already in mm ---------")
	
                #Lead the Vtu file from mesh-complete folder
		MeshVTK=ReadVTUFile(self.Args.InputFolder+"/mesh-complete.mesh.vtu")
 
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


		self.IDs["inlet_id"]=[1]
		self.IDs["outlet_ids"]=[]
		OutletSurfaceArea=[]
		for BoundaryFile_ in BoundaryFiles:
			print ("------ Looping over: %s"%BoundaryFile_)
			
			BoundarySurface_=ReadVTPFile(BoundaryFile_)
			BoundaryNodeIds+=[BoundarySurface_.GetPointData().GetArray("GlobalNodeID").GetValue(i) for i in range(BoundarySurface_.GetNumberOfPoints())]
			BoundaryEntityIds+=[counter]*BoundarySurface_.GetNumberOfPoints()

			#Add the centroid and surface area of the boundary file to the self.IDs dictionary
			SurfaceArea_=GetSurfaceArea(BoundarySurface_)
			Centroid_=GetCentroid(BoundarySurface_)
			if counter>=2:
				self.IDs[BoundaryFile_.split("/")[-1].replace(".vtp","")]=Centroid_
				self.IDs[BoundaryFile_.split("/")[-1].replace(".vtp","_area")]=SurfaceArea_
			
			#All of the outlets
			if counter>2: 
				self.IDs["outlet_ids"].append(counter-1)
				OutletSurfaceArea.append(SurfaceArea_)				
			counter+=1
			
		self.IDs["mean_flow_rate"]=1.00

		AreaRatio=np.array(OutletSurfaceArea)/np.sum(OutletSurfaceArea)
		self.IDs["area_ratio"]=np.ndarray.tolist(AreaRatio)

		print ("--- Writing the Boundary IDS to: %s"%self.Args.OutputFileName+"_info.json")
		with open(self.Args.OutputFileName+"_info.json", 'w') as convert_file:
			convert_file.write(json.dumps(self.IDs))


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
		print ("--- Writing the file in VTU format: %s"%self.Args.OutputFileName+".vtu")
		WriteVTUFile(self.Args.OutputFileName+".vtu",MeshVTK)

		#print ("--- Scaling the mesh from cm to mm since SimVascular is usually in cm")
		#print ("--- Using the Scaling Factor: %.02f"%self.Args.ScalingFactor)
		#os.system("vmtkmeshscaling -ifile %s -ofile %s -scale %f"%(self.Args.OutputFileName+".vtu",self.Args.OutputFileName+".vtu",self.Args.ScalingFactor))

		#Write the Mesh in xml format
		print ("--- Writing the file in XML format: %s"%self.Args.InputFolder+".xml")
		os.system("vmtkmeshwriter -ifile %s.vtu -ofile %s.xml -entityidsarray CellEntityIds"%(self.Args.OutputFileName,self.Args.OutputFileName))

	def XMLMeshWriter(self,Mesh,OutputFileName):
		MeshWriter=vmtkscripts.vmtkMeshWriter()
		MeshWriter.Mesh=Mesh
		MeshWriter.Format="dolfin"
		MeshWriter.GuessFormat=0
		MeshWriter.OutputFileName=OutputFileName
		MeshWriter.Compressed=1
		MeshWriter.Execute()



if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take a mesh-complete folder from SimVascular and write a dolfin mesh file. The CellEntityIds 0 is for the volume (Tetrahedron), 1 for the mesh wall, 2 for inlet, and 3....N for outlets.")
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The path to mesh-complete folder from SimVascular")
<<<<<<< HEAD
	parser.add_argument('-ScalingFactor', '--ScalingFactor', type=int, required=False,help="Scale the mesh from cm to mm.")

	#Output Filename 
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName",help="The output file in which to store the dolfin mesh.")

	args=parser.parse_args()
	OasisMeshWriterForSimVascular(args).Main()



 
