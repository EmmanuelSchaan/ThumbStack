#!/bin/bash
#SBATCH -J tszksz_pact
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use haswell/knl nodes
#SBATCH -t 11:59:59  #00:29:00  #09:59:59  #30:00:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eschaan@lbl.gov
#SBATCH -o /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/tszksz_pact.out
#SBATCH -e /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/log/tszksz_pact.err


cd /global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/
source ~/python_profile.sh


python driver_tszksz_pactf150night_2019_03_11.py
python driver_tszksz_pactf090night_2019_03_11.py

#python driver_tszksz_tilec_pact_cmbksz_d56.py
#python driver_tszksz_tilec_pact_cmbksz_boss.py
#python driver_tszksz_tilec_pact_comptony_boss.py
#python driver_tszksz_tilec_pact_comptonydeprojcib_boss.py

#python driver_tszksz_planck_smicacmb.py
#python driver_tszksz_planck_smicacmbnosz.py
#
#python driver_tszksz_planck545.py

