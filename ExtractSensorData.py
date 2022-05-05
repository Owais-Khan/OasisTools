import os, sys
import vtk, numpy, h5py
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
import argparse
from glob import glob
from utilities import *

class ExtractSensorData():
	def __init__(self,Args):
		self.Args=Args
	
		if self.Args.OutFileName is None:
			self.Args.OutFileName=self.Args.VelocityFileName.replace(".vtu","_%d_sensors.txt"%self.Args.Sensors)	

	def Main(self):
		#Read the input velocoity
		print ("--- Loading Velocity: %s"%self.Args.VelocityFileName)
		Velocity=ReadVTUFile(self.Args.VelocityFileName)

		#Extract the Surface
		print ("--- Extracting the Surface")
		Surface=ExtractSurface(Velocity)
		
		#Get SurfaceCoordinates
		SurfaceCoords=np.array([Surface.GetPoint(i) for i in range(Surface.GetNumberOfPoints())])

		#Get Number of Surface Points
		Npts_v=Velocity.GetNumberOfPoints()
		Npts_s=len(SurfaceCoords)
		Incr_=int((Npts_v-Npts_s)/self.Args.Sensors)

		print ("--- Number of volume points : %d"%Npts_v)

		#Loop over the volume points to get the number of sensors
		print ("--- Extract the Sensor Points")
		Coords=np.zeros(shape=(self.Args.Sensors,3))
		VelocityArray=np.zeros(shape=(self.Args.Sensors,3))
		PressureArray=np.zeros(self.Args.Sensors)
		counter=0
		for i in range(Npts_s,Npts_v,Incr_):
			if counter==self.Args.Sensors: break
			Coords[counter,:]=Velocity.GetPoint(i)	
			VelocityArray[counter,0]=Velocity.GetPointData().GetArray("Velocity").GetValue(i*3)
			VelocityArray[counter,1]=Velocity.GetPointData().GetArray("Velocity").GetValue(i*3+1)
			VelocityArray[counter,2]=Velocity.GetPointData().GetArray("Velocity").GetValue(i*3+2)
			counter+=1
	
		counter=0
		for i in range(Npts_s,Npts_v,Incr_):
			if counter==self.Args.Sensors: break
			PressureArray[counter]  =Velocity.GetPointData().GetArray("Pressure").GetValue(i)
			counter+=1
		
		#Write the sensor data
		print ("--- Saving the Outfile: %s"%self.Args.OutFileName)
		self.WriteSensors(VelocityArray,PressureArray,Coords,self.Args.OutFileName)

	def WriteSensors(self,Velocity,Pressure,Coords,OutFileName):
		#Record the sensor data
		outfile=open(OutFileName,'w')
		outfile.write("X Y Z Velocity0 Velocity1 Velocity2 Pressure\n")
		for i in range(self.Args.Sensors):
			outfile.write("%.08f %.08f %.08f %.08f %.08f %.08f %.08f\n"%(Coords[i,0],Coords[i,1],Coords[i,2],Velocity[i,0],Velocity[i,1],Velocity[i,2],Pressure[i]))
		outfile.close()	
			
		
			

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will extract the sensor data (Velocity, Pressure) for PINNS calculations. You need to provide the # of sensors to extract.")
        
	parser.add_argument('-VelocityFileName', '--VelocityFileName', type=str, required=False, dest="VelocityFileName", help="The file that contains velocity and pressure in vtu format")
	
	parser.add_argument('-Sensors', '--Sensors', type=int, required=True, dest="Sensors", help="The number of sensors to gather from the data.")
        
	parser.add_argument('-OutFileName', '--OutFileName', type=str, required=False, dest="OutFileName", help="The filename to store the output sensor data.")

        
	args=parser.parse_args()
        
	ExtractSensorData(args).Main()


