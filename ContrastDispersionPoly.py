import vtk
import numpy as np
import os
import glob
import argparse
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
from utilities import ReadXDMFFile, vtk_to_numpy
class ContrastDispersionPoly():
    def __init__(self,args) -> None:
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
        plane.SetOrigin(min_val, 0., 0.)
        slicer = vtk.vtkExtractGeometry()
        slicer.SetExtractInside(1)
        slicer.ExtractBoundaryCellsOn()
        slicer.SetImplicitFunction(plane)
        time_array = []
        lag = 2
        time_steps = np.arange(0,N,int(N/NFiles))
        for i in time_steps[lag:]:
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
        t = time_steps[lag:]/1000
        t = t.reshape(-1,1)
        x = np.arange(min_val, max_val, interval).reshape(-1,1)
        time_array = np.array(time_array).reshape(-1,1)
        space_array = np.array(space_array).reshape(-1,1)
        poly = PolynomialFeatures(degree = 1, include_bias=False)
        poly_features = poly.fit_transform(t)
        model = LinearRegression()
        model.fit(poly_features,time_array)
        pred = model.predict(poly_features)
        plt.scatter(t, time_array)
        plt.plot(t, pred, color = 'red')
        plt.show()
        dc_dt = model.coef_[0][0]
		# fitting a 2nd order polynomial on spatial contrast array
        poly = PolynomialFeatures(degree=2, include_bias=False)
        poly_features = poly.fit_transform(x)
        model.fit(poly_features, space_array)
        coef2 = model.coef_[0][1]
        coef1 = model.coef_[0][0]
        d2c_dx2 = 2*coef2
        print(d2c_dx2)
        dc_dx = np.array(coef2*x + coef1)
        dc_dx = dc_dx[0][0]
        print(dc_dx)
        pred = model.predict(poly_features)
        plt.scatter(x, space_array)
        plt.plot(x, pred, color = 'red')
        plt.show()
        D = 0.04
        velocity = (D*d2c_dx2 - dc_dt)/dc_dx
        print('The fluid velocity is: ',abs(velocity))
if __name__=='__main__':
    parser = argparse.ArgumentParser(description="This script implement the contrast dispersion pipeline")
    parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest='InputFolder')
    args = parser.parse_args()
    ContrastDispersionPoly(args).main()