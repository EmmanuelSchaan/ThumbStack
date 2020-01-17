#!/bin/bash
#SBATCH -J gengrf
#SBATCH -N 1
#SBATCH -q debug  #regular
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use haswell/knl nodes
#SBATCH -t 00:29:00  #47:59:59 #05:00:00  #27:00:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eschaan@lbl.gov
#SBATCH -o /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/gengrf_mocks.out
#SBATCH -e /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/gengrf_mocks.err


cd /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/
source ~/python_profile.sh


# Compute everything in chunks
python generate_mocks_grf_f150_night_2019_03_11.py 0 50
#python generate_mocks_grf_f150_night_2019_03_11.py 50 50
#python generate_mocks_grf_f150_night_2019_03_11s.py 100 50
#python generate_mocks_grf_f150_night_2019_03_11s.py 150 50
#python generate_mocks_grf_f150_night_2019_03_11s.py 200 50
#python generate_mocks_grf_f150_night_2019_03_11s.py 250 50
#python generate_mocks_grf_f150_night_2019_03_11s.py 300 50
#python generate_mocks_grf_f150_night_2019_03_11s.py 350 50
#python generate_mocks_grf_f150_night_2019_03_11s.py 400 50

# Run everything together (after modifying the code not to recompute)
#python generate_mocks_grf_f150_night_2019_03_11s.py 0 400
