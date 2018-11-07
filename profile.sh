## To install pixell
## git clone the repository
## First, only load the following:
#module load python/2.7-anaconda-4.4
#module load intel
#module unload openmpi
#module unload impi
#module load gcc/4.9.3
## Then, go inside the pixell folder and run:
#LD=icc CC=icc FC=ifort F77=ifort python setup.py build_ext -i --fcompiler=intelem --compiler=intelem
## After this, you can load the other modules, and activate a conda environment.
## The reason for not activating a conda environment right away is
## that pip install then tries to use the conda fortran compiler,
## as opposed to the module loaded gcc

# Compilers
module load intel
module load gcc/4.9.3

# OpenMPI
module load openmpi

# FFTW
if [ "$NERSC_HOST" = "cori" ]; then
   module load cray-fftw/3.3.6.2
   #module load fftw/3.3.4.6
else
   module load cray-fftw/3.3.6.3
   #module load fftw/3.3.4.10
fi

# python
module load python/2.7-anaconda-4.4
export PATH=$PYTHONUSERBASE/bin:$PATH

# conda environment for python
source activate pixell

# Add path to pixell code to the python path
export PYTHONPATH=$PYTHONPATH:/global/homes/e/eschaan/local/pixell
