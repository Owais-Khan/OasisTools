#This program is written by Owais Khan
#This program generates the 7 points for LCCA

import numpy
import os, sys
import math
from scipy import interpolate
from vmtk import vtkvmtk
import vtk
import argparse

class OasisCreateSurface():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			FileName_=self.Args.InputFileName.split("/")[-1]
			self.Args.OutputFolder=self.Args.InputFileName.replace(FileName_,"")+"Results_Surface"
			os.system("mkdir %s"%self.Args.OutputFolder)

	def Main(self):
		############################################################
		############### Generate Surface ###########################
		############################################################
		CentrePointList=open(self.Args.InputFileName,'r')

		CoarseCentreWriter = vtk.vtkXMLPolyDataWriter()
		CoarseCentreWriter.SetFileName('%s/step1_Coarse_CL.vtp'%self.Args.OutputFolder)

		FineCentreWriter = vtk.vtkXMLPolyDataWriter()
		FineCentreWriter.SetFileName('%s/step2_Fine_CL.vtp'%self.Args.OutputFolder)

		ObjectImageWriter = vtk.vtkXMLImageDataWriter()
		ObjectImageWriter.SetFileName('%s/step3_image.vti'%self.Args.OutputFolder)

		CentrePoints=vtk.vtkPoints()
		CentreCellIDs=vtk.vtkIdList()
		CentreCellsArray=vtk.vtkCellArray()
		RadsVTK = vtk.vtkDoubleArray()
		RadsVTK.SetName("MaximumInscribedSphereRadius")
		n=0
		for FileLine in CentrePointList:
			#print [float(FileLine.split()[0]), float(FileLine.split()[1]), float(FileLine.split()[2])]
        		
			CentreCellIDs.InsertNextId(CentrePoints.InsertNextPoint([float(FileLine.split()[0]), float(FileLine.split()[1]), float(FileLine.split()[2])]))
     
			RadsVTK.InsertNextValue(float(FileLine.split()[3]))
			n+=1



		CentreCellsArray.InsertNextCell(CentreCellIDs)
		CentreLineData=vtk.vtkPolyData()
		CentreLineData.SetPoints(CentrePoints)
		CentreLineData.SetLines(CentreCellsArray)
		CentreLineData.GetPointData().AddArray(RadsVTK)

		CoarseCentreWriter.SetInputData(CentreLineData)
		CoarseCentreWriter.Write()

		x=numpy.zeros(n,float)
		y=numpy.zeros(n,float)
		z=numpy.zeros(n,float)
		rad=numpy.zeros(n,float)
		temp=[0.0,0.0,0.0]
		for i in range (0,n):
			CentrePoints.GetPoint(i, temp)
			x[i]=temp[0]
			y[i]=temp[1]
			z[i]=temp[2]
			rad[i]=RadsVTK.GetValue(i)
		
		tck,u = interpolate.splprep([x,y,z,rad],k=1, task=0, s=0)
		unew = numpy.arange(0,1.005,0.001)
		out = interpolate.splev(unew,tck)

		CentrePointsHigh=vtk.vtkPoints()
		CentreCellIDsHigh=vtk.vtkIdList()
		CentreCellsArrayHigh=vtk.vtkCellArray()
		RadsHighVTK = vtk.vtkDoubleArray()
		RadsHighVTK.SetName("MaximumInscribedSphereRadius")

		for i in range(0,unew.size):
			CentreCellIDsHigh.InsertNextId(CentrePointsHigh.InsertNextPoint(out[0][i],out[1][i], out[2][i]))
			RadsHighVTK.InsertNextValue(out[3][i])

		CentreCellsArrayHigh.InsertNextCell(CentreCellIDsHigh)
		CentreLineDataHigh=vtk.vtkPolyData()
		CentreLineDataHigh.SetPoints(CentrePointsHigh)
		CentreLineDataHigh.SetLines(CentreCellsArrayHigh)
		CentreLineDataHigh.GetPointData().AddArray(RadsHighVTK)

		FineCentreWriter.SetInputData(CentreLineDataHigh)
		FineCentreWriter.Write()

		PBModeller = vtkvmtk.vtkvmtkPolyBallModeller()
		PBModeller.SetInputData(CentreLineDataHigh)
		PBModeller.SetRadiusArrayName("MaximumInscribedSphereRadius")
		PBModeller.UsePolyBallLineOff()
	
		#Find the bounding box distance
		DiagDist=numpy.sqrt( (max(x)-min(x))**2 + (max(y)-min(y))**2 +(max(z)-min(z))**2)
		Resolution=DiagDist/self.Args.SampleDimensions
	
		#Define and isotropic resolution in X, Y and Z directions
		Nx=int((max(x)-min(x))/Resolution)
		Ny=int((max(y)-min(y))/Resolution)
		Nz=int((max(z)-min(z))/Resolution)
		
		if Nx==0: Nx=max(Nx,Ny,Nz)		
		if Ny==0: Ny=max(Nx,Ny,Nz)		
		if Nz==0: Nz=max(Nx,Ny,Nz)		

		PBModeller.SetSampleDimensions([Nx,Ny,Nz])
		PBModeller.Update()

		ObjectImageWriter.SetInputData(PBModeller.GetOutput())
		ObjectImageWriter.Write()
		ObjectImageWriter.Update()
		os.system('vmtkmarchingcubes -ifile %s/step3_image.vti -ofile %s/step4_surface.vtp'%(self.Args.OutputFolder,self.Args.OutputFolder))
		os.system('vmtksurfaceviewer -ifile %s/step4_surface.vtp'%self.Args.OutputFolder)

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script can generate surfaces with tubular geometries based on coordiantes provided in a text file")
        
	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The filename where X Y Z Radius data is stored")
	
	parser.add_argument('-SampleDimensions', '--SampleDimensions', type=str, required=False, default=100, dest="SampleDimensions",help="The isotropic dimens isotropic dimensions in x,y and z directions. Default is 100")
        
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder where the results will be stored")

        
	args=parser.parse_args()
	OasisCreateSurface(args).Main()


