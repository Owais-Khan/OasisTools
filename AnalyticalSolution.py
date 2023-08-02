#This script was written by Anahita Seresti on July 28, 2023

from dolfin import *
import argparse

class AnalyticalSolution():
	def __init__(self,args):
		self.Args = args
	def main(self):
		Mesh = IntervalMesh(self.Args.NofElements, self.Args.minBound, self.Args.maxBound)
		c1 = Constant(self.Args.ConstantVal1)
		c2 = Constant(self.Args.ConstantVal2)
		V = Constant(self.Args.MeanVelocity)
		D = Constant(self.Args.DiffusionCoef)
		e = Constant(exp(V/D))
		U = FunctionSpace(Mesh,"Lagrange",1)
		c = interpolate(Expression("(c2-c1)/(e-1)*exp(V/D*x[0])+(c1*e-c2)/(e-1)",c1 = c1, c2 = c2, e = e, V = V, D = D, degree = 2),U)
		ofile = XDMFFile("./Concentration1D.xdmf")
		ofile.write(c)
		return c


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "This script calculates concentration based on the 1D analytical solution for the steady advection diffusion")
	parser.add_argument('-NofElements', '--NofElements', type=int, required=False, default = 100, dest = "NofElements")
	parser.add_argument('-ConstantVal1', '--ConstantVal1', type=float, required=False, default = 1.0, dest = "ConstantVal1")
	parser.add_argument('-ConstantVal2', '--ConstantVal2', type=float, required=False, default = 0.1, dest = "ConstantVal2")
	parser.add_argument('-minBound', '--minBound', type=float, required=False, default = 0., dest = "minBound")
	parser.add_argument('-maxBound', '--maxBound', type=float, required=False, default = 10., dest = "maxBound")
	parser.add_argument('-MeanVelocity', '--MeanVelocity', type=float, required=False, default = 50.0, dest = "MeanVelocity")
	parser.add_argument('-DiffusionCoef', '--DiffusionCoef', type=float, required=False, default = 1.0, dest = "DiffusionCoef")

	args = parser.parse_args()
	AnalyticalSolution(args).main()
