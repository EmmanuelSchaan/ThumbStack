#!/bin/bash
#SBATCH -J testendtoend
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use haswell/knl nodes
#SBATCH -t 04:00:00  #30:00:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eschaan@lbl.gov
#SBATCH -o /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/test_endtoend.out
#SBATCH -e /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/test_endtoend.err


cd /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/
source ~/python_profile.sh


# Old CMB maps
#python generate_mocks_grf_f150_daynight_2018_08_10.py

# New CMB maps
python driver_test_endtoend_carmanu_tszksz_f150night_2019_03_11.py

