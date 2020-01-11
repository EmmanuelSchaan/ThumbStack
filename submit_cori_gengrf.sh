#!/bin/bash
#SBATCH -J gengrf
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use haswell/knl nodes
#SBATCH -t 20:00:00 #05:00:00  #27:00:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eschaan@lbl.gov
#SBATCH -o /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/gengrf_mocks.out
#SBATCH -e /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/gengrf_mocks.err


cd /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/
source ~/python_profile.sh


# Old CMB maps
#python generate_mocks_grf_f150_daynight_2018_08_10.py

# New CMB maps
#python generate_mocks_grf_f150_night_2019_03_11.py 101 2

# Compute everything in chunks
#python generate_mocks_grf_f150_night_2019_03_11.py 100 100
#python generate_mocks_grf_f150_night_2019_03_11.py 200 100
#python generate_mocks_grf_f150_night_2019_03_11.py 300 100

# Run everything together (after modifying the code not to recompute)
python generate_mocks_grf_f150_night_2019_03_11.py 0 400
