import numpy as np
import vtk
from utilities import vtk_to_numpy, ReadXDMFFile
import os
import glob
import argparse
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
class ContrastDispersion():
	def __init__(self,args):
		self.Args = args
	def main(self):
		# Get the input Folder and take 10 timesteps with constant intervals turn them to vtu file and read them
		filenames = glob.glob(f'{self.Args.InputFolder}/*.xdmf')
		filenames = sorted(filenames)
		N = len(filenames)
		NFiles = 10
		Mesh = {i:ReadXDMFFile(filenames[i]) for i in range(0,N,int(N/NFiles))}
		# slice: take the first slice of all time steps
		min_val = Mesh[0].GetBounds()[0]
		max_val = Mesh[0].GetBounds()[1]
		normal = (1., 0., 0.)
		plane = vtk.vtkPlane()
		plane.SetNormal(normal)
		plane.SetOrigin(Mesh[0].GetBounds()[0], 0., 0.)
		slicer = vtk.vtkExtractGeometry()
		slicer.SetExtractInside(1)
		slicer.ExtractBoundaryCellsOn()
		slicer.SetImplicitFunction(plane)
		time_array = []
		for i in range(0,N,int(N/NFiles)):
			slicer.SetInputData(Mesh[i])
			slicer.Update()
			array_ = slicer.GetOutput().GetPointData().GetArray('f_22')
			array_ = vtk_to_numpy(array_)
			time_array.append(np.average(array_))
			peak_val = i
		# take the slices over the last mesh
		NSlices = 100
		interval = (max_val - min_val)/NSlices
		slicer.SetInputData(Mesh[peak_val])
		space_array = []
		for i in range(NSlices):
			point = min_val + i*interval
			plane.SetOrigin(point, 0., 0.)
			slicer.SetImplicitFunction(plane)
			slicer.Update()
			array = slicer.GetOutput().GetPointData().GetArray('f_22')
			array = vtk_to_numpy(array)
			space_array.append(np.average(array))
		# fit line on time and space domain
		t = np.arange(0,N,int(N/NFiles))/1000
		x = np.arange(min_val, max_val, interval)
		time_array = np.array(time_array)
		space_array = np.array(space_array)
		model = LinearRegression()
		model.fit(t.reshape(-1,1),time_array.reshape(-1,1))
		pred = model.predict(t.reshape(-1,1))
		plt.scatter(t.reshape(-1,1), time_array.reshape(-1,1))
		plt.plot(t.reshape(-1,1), pred.reshape(-1,1), color = 'red')
		plt.show()
		dc_dt = model.coef_[0][0]
		model.fit(x.reshape(-1,1),np.array(space_array).reshape(-1,1))
		dc_dx = model.coef_[0][0]
		pred = model.predict(x.reshape(-1,1))
		plt.scatter(x.reshape(-1,1), space_array.reshape(-1,1))
		plt.plot(x.reshape(-1,1), pred.reshape(-1,1), color = 'red')
		plt.show()
		velocity = dc_dt/dc_dx
		print(dc_dt,dc_dx,abs(velocity))

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="This script implement the contrast dispersion pipeline")
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest='InputFolder')
	args = parser.parse_args()
	ContrastDispersion(args).main()
