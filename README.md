# OasisTools
The scripts in this repository provide useful tools to facilitate CFD simulations using the Oasis solver (https://github.com/mikaem/Oasis). You will need the following libraries to use these tools:
1. vtk
2. Fenics

If you have a mac, I recommend installing VTK first and then Fenics using the following commands:

```console
foo@bar:~$ conda create -n Fenics
foo@bar:~$ conda activate Fenics
foo@bar:~$ conda install vtk
foo@bar:~$ conda install -c conda-forge fenics
```

To get help for any script, please type:
```console
foo@bar:~$ python [ScriptName.py] -h
```
## Converting SimVascular Mesh to Oasis readable Mesh
SimVascular software (see simvascular.github.io) can be used to segment and generate volumetric meshes. SimVascular produces a folder called "mesh-complete", which contains the volumetric, wall and cap meshes. Prior to generating the mesh in SimVascular, please ensure that the inflow cap is labelled as "inflow". 

The following script can be used to produce a VTU file that can be visualized in Paraview and .xml.gz file that can be read into Oasis CFD solver. Here is how to use the script.

```console
foo@bar:~$ python OasisMeshWriterForSimVascular.py -InputFolder /path/to/mesh-complete/ 
```
You may define the following parameters:
The script will output two files in the same folder as the /path/to/ (i.e. where mesh-complete folder is located). 
1. mesh-complete.vtu: A volumetric mesh file that can be read into paraview along with all of the boundary ids.
2. mesh-complete.xml.gz: A volumetric msh file that can be read into Oasis. Please refer to VaMPY documentation (https://github.com/KVSlab/VaMPy) on how to conduct simulations with the Oasis CFD solver.

**Note**: The CellEntityIds tag in the mesh correspond to: 0=volumetric mesh, 1=wall mesh, 2=inflow and 3...N for outflows.

**Note**: The VTU file is converted in .xml.gz file using the vmtkMeshWriter script. However, the conda installing of vmtk will give you an error. You will need to update the vmtkMeshWriter script manually to corrct this. Please open the vmtkmeshwriter.py file, go to line 264, and change ```file = open(self.OutputFileName,'r')``` to ```file = open(self.OutputFileName,'rb')```. If you install vmtk using conda, the script will be located at: ```/Users/[USERNAME]/miniconda3/envs/vmtk/lib/python3.6/site-packages/vmtk/vmtkmeshwriter.py```


 


