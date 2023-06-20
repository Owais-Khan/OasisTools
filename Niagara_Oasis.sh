#!/bin/bash

# Name of your job
#SBATCH --job-name=Oasis
#SBATCH --partition=compute

# Specify the name of the output file. The %j specifies the job ID
#SBATCH --output=Oasis.o%j

# Specify the name of the error file. The %j specifies the job ID
#SBATCH --error=Oasis.e%j

# The walltime you require for your job
#SBATCH --time=01:00:00

# Job priority. Leave as normal for now
#SBATCH --qos=normal

# Number of nodes are you requesting for your job. You can have 40 processors per node
#SBATCH --nodes=1

# Number of processors per node
#SBATCH --ntasks-per-node=40

# Send an email to this address when your job starts and finishes
#SBATCH --mail-user=owaiskhan@ryerson.ca
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Name of the executable you want to run on the cluster
module purge; module load intel/2019u3 intelmpi/2019u3 boost/1.69.0 hdf5 netcdf/4.6.3 petsc/3.10.5 trilinos/12.12.1 intelpython3/2019u3 fenics/2019.1.0; module load hdf5

mpirun -np 40 /home/k/khanmu11/khanmu11/Softwares/Oasis/bin/oasis NSfracStep problem=Artery mesh_path=./mesh-complete.xml.gz 

