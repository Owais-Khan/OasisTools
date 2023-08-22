import vtk
import numpy as np
import glob
import argparse
from matplotlib.pyplot import plot as plt
from utilities import ReadXDMFFile, vtk_to_numpy
class VelocityPointByPoint():
    def __init__(self, args) -> None:
        self.Args = args
    def main(self):
        pass
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script takes the folder at the output of the OasisAdvectionDiffusion.py")
    parser.add_argument('-InputFolder', '--InputFolder', required=True, type=str, dest='InputFolder', help="Input Folder containing the Concentration Files in .XDMF format")
