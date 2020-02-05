#!/bin/bash
#SBATCH -J gengrf
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use haswell/knl nodes
#SBATCH -t 47:59:59  #47:59:59 #05:00:00  #27:00:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eschaan@lbl.gov
#SBATCH -o /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/gengrf_mocks.out
#SBATCH -e /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/gengrf_mocks.err


cd /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/
source ~/python_profile.sh


# Compute everything in chunks
#python generate_mocks_grf_f150_night_2019_03_11.py 0 50
#python generate_mocks_grf_f150_night_2019_03_11.py 50 50
#python generate_mocks_grf_f150_night_2019_03_11.py 100 50
#python generate_mocks_grf_f150_night_2019_03_11.py 150 50
#python generate_mocks_grf_f150_night_2019_03_11.py 200 50
#python generate_mocks_grf_f150_night_2019_03_11.py 250 50
#python generate_mocks_grf_f150_night_2019_03_11.py 300 50
#python generate_mocks_grf_f150_night_2019_03_11.py 350 50
#python generate_mocks_grf_f150_night_2019_03_11.py 400 50
#
#python generate_mocks_grf_f150_night_2019_03_11.py 494 50
#python generate_mocks_grf_f150_night_2019_03_11.py 500 50
#python generate_mocks_grf_f150_night_2019_03_11.py 550 50
#python generate_mocks_grf_f150_night_2019_03_11.py 600 50
#python generate_mocks_grf_f150_night_2019_03_11.py 650 50
#python generate_mocks_grf_f150_night_2019_03_11.py 700 50
#python generate_mocks_grf_f150_night_2019_03_11.py 750 50

# Run everything together (after modifying the code not to recompute)
#python generate_mocks_grf_f150_night_2019_03_11.py 0 400
python generate_mocks_grf_f150_night_2019_03_11.py 0 800



# just to fix some of the mocks
#python generate_mocks_grf_f150_night_2019_03_11.py 652 149

