# This script was developed by Anahita Seresti on July 28-Aug 1
import argparse
import numpy as np
import vtk
import os
import sys
from utilities import ReadVTUFile, vtk_to_numpy
from OasisAnalyticalSolution import AnalyticalSolution
class SimulationVsAnalytical():
	def __init__(self,args):
		self.Args = args
	def main(self):
		NSlice = 100
		Mesh3D = ReadVTUFile(self.Args.InputMeshFile)
		min_val = Mesh3D.GetBounds()[0]
		max_val = 1.#Mesh3D.GetBounds()[1]
		interval = (max_val - min_val)/NSlice
		normal = (1.,0.,0.)
		plane = vtk.vtkPlane()
		plane.SetNormal(normal)
		slicer = vtk.vtkExtractGeometry()
		slicer.SetInputData(Mesh3D)
		slicer.SetExtractInside(1)
		slicer.ExtractBoundaryCellsOn()
		plane.SetOrigin(min_val, 0., 0.)
		slicer.SetImplicitFunction(plane)
		slicer.Update()
		array = slicer.GetOutput().GetPointData().GetArray('f_22')
		array = vtk_to_numpy(array)
		ConstantVal2 = np.average(array)
		plane.SetOrigin(min_val+NSlice*interval, 0., 0.)
		slicer.SetImplicitFunction(plane)
		slicer.Update()
		array = slicer.GetOutput().GetPointData().GetArray('f_22')      
		array = vtk_to_numpy(array)
		ConstantVal1 = np.average(array)
		self.Args.ConstantVal1 = ConstantVal1; self.Args.ConstantVal2 = ConstantVal2; self.Args.NofElements = NSlice; self.Args.minBound = 0.; self.Args.maxBound = 1.#max_val
		Mesh1D = AnalyticalSolution(self.Args)
		c = Mesh1D.main()
		analytical_array = c.vector().get_local()
		mean_error = 0
		ofile = './validation.txt'
		with open(ofile, 'w') as file:
			file.write('Simulation, Analytical, Error(%) \n')
			for i in range(NSlice+1):
				point = min_val + i*interval
				print(f'Slicing the mesh on x = {point}')
				plane.SetOrigin(point,0.,0.)
				slicer.SetImplicitFunction(plane)
				slicer.Update()
				array = slicer.GetOutput().GetPointData().GetArray('f_22')
				array = vtk_to_numpy(array)
				error = abs(analytical_array[i]-np.average(array))/analytical_array[i]*100
				mean_error += error**2
				file.write(f'{np.average(array)}, {analytical_array[i]}, {error}\n')
			file.write(f'Mean Squared Error = {round((mean_error/NSlice)**0.5,3)}%')
		file.close()
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description = "This script takes a mesh file and outputs the average value of the plane cuts along the Mesh")
	parser.add_argument('-InputMeshFile', '--InputMeshFile',type=str, required=True, dest='InputMeshFile')
	parser.add_argument('-NofElements', '--NofElements', type=int, required=False, dest = "NofElements")
	parser.add_argument('-ConstantVal1', '--ConstantVal1', type=float, required=False, dest = "ConstantVal1")
	parser.add_argument('-ConstantVal2', '--ConstantVal2', type=float, required=False, dest = "ConstantVal2")
	parser.add_argument('-minBound', '--minBound', type=float, required=False, dest = "minBound")
	parser.add_argument('-maxBound', '--maxBound', type=float, required=False, dest = "maxBound")
	parser.add_argument('-MeanVelocity', '--MeanVelocity', type=float, required=False, default = 50.0, dest = "MeanVelocity")
	parser.add_argument('-DiffusionCoef', '--DiffusionCoef', type=float, required=False, default = 1.0, dest = "DiffusionCoef")
	args = parser.parse_args()
	SimulationVsAnalytical(args).main()
