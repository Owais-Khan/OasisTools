import vtk
import numpy as np
import glob
import argparse
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
from utilities import ReadXDMFFile, vtk_to_numpy
class VelocityPointByPoint():
    def __init__(self, args) -> None:
        self.Args = args
        filenames = glob.glob(f'{self.Args.InputFolder}/*.xdmf')
        filenames = sorted(filenames)
        self.N = len(filenames)
        self.Mesh = {i:ReadXDMFFile(filenames[i]) for i in range(0,self.N,int(self.N/self.Args.NFiles))}
		# slice: take the first slice of all time steps
        self.min_val = self.Mesh[0].GetBounds()[0]
        self.max_val = self.Mesh[0].GetBounds()[1]
        self.interval = (self.max_val - self.min_val)/self.Args.NSlices
    def main(self):
        normal = (1., 0., 0.)
        plane = vtk.vtkPlane()
        plane.SetNormal(normal)
        plane.SetOrigin(self.min_val, 0., 0.)
        slicer = vtk.vtkExtractGeometry()
        slicer.SetExtractInside(1)
        slicer.ExtractBoundaryCellsOn()
        slicer.SetImplicitFunction(plane)
        time_array = []
        lag = 0
        time_steps = np.arange(0,self.N,int(self.N/self.Args.NFiles))
        for i in time_steps[lag:]:
            slicer.SetInputData(self.Mesh[i])
            slicer.Update()
            array_ = slicer.GetOutput().GetPointData().GetArray('f_22')
            array_ = vtk_to_numpy(array_)
            time_array.append(np.average(array_))
            peak_val = i
        # fit a linear model on the upslope
        t = time_steps[lag:]/1000
        t = t.reshape(-1,1)
        time_array = np.array(time_array).reshape(-1,1)
        poly = PolynomialFeatures(degree = 1, include_bias=False)
        poly_features = poly.fit_transform(t)
        model = LinearRegression()
        model.fit(poly_features,time_array)
        pred = model.predict(poly_features)
        plt.scatter(t, time_array)
        plt.plot(t, pred, color = 'red')
        plt.show()
        dc_dt = model.coef_[0][0]
		# take the slices over the last mesh
        self.interval = (self.max_val - self.min_val)/self.Args.NSlices
        slicer.SetInputData(self.Mesh[peak_val])
        space_array = []
        velocity = []
        for i in range(self.Args.NSlices):
            point = self.min_val + i*self.interval
            plane.SetOrigin(point, 0., 0.)
            slicer.SetImplicitFunction(plane)
            slicer.Update()
            array = slicer.GetOutput().GetPointData().GetArray('f_22')
            array = vtk_to_numpy(array)
            space_array.append(np.average(array))
        dc_dx, d2c_dx2 = VelocityPointByPoint(args).FitSpatialModel(space_array)
        for i in range(self.Args.NSlices):
            velocity.append((self.Args.D*d2c_dx2 - dc_dt)/dc_dx[i])
        print(velocity)
        print('The fluid velocity is: ',abs(np.mean(np.array(velocity))))
    def FitSpatialModel(self,space_array):
        x = np.arange(self.min_val, self.max_val, self.interval).reshape(-1,1)
        space_array = np.array(space_array).reshape(-1,1)
        poly = PolynomialFeatures(degree=2, include_bias=False)
        poly_features = poly.fit_transform(x)
        model = LinearRegression()
        model.fit(poly_features, space_array)
        coef2 = model.coef_[0][1]
        coef1 = model.coef_[0][0]
        d2c_dx2 = 2*coef2
        #print(d2c_dx2)
        dc_dx = np.array(coef2*x + coef1)
        #dc_dx = dc_dx[0][0]
        return dc_dx, d2c_dx2
        #print(dc_dx)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script takes the folder at the output of the OasisAdvectionDiffusion.py")
    parser.add_argument('-InputFolder', '--InputFolder', required=True, type=str, dest='InputFolder', help="Input Folder containing the Concentration Files in .XDMF format")
    parser.add_argument('-D', '--D', required=False, default=0.04, dest='D', type= float, help='The diffusion coefficient')
    parser.add_argument('-NFiles', '--NFiles', required=False, default=10, type=int, dest='NFiles', help="number of temporal points")
    parser.add_argument('-NSlices', '--NSlices', required=False, default=100, type=int, dest='NSlices', help="The number of the slices along the vessels")
    args = parser.parse_args()
    Object = VelocityPointByPoint(args)
    Object.main()