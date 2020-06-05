
# create sub-folders if needed
#mkdir cmass
#mkdir lowz
#mkdir pipeline

####################################################################################
# footprint of BOSS vs ACT
#! Redo for Kendrick
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/catalog/cmass_mariana/footprint_cmass_mariana.pdf ./

# summary dn/dz
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/catalog/cmass_mariana/summary_dndz.pdf ./


####################################################################################
# CMASS

# summary ksz
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/cmass_kendrick_pactf150daynight20200228maskgal60/summary_ksz_150_90_cmass_kendrick.pdf ./cmass/
# joint cov 150-90
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/cmass_kendrick_pactf150daynight20200228maskgal60/cor_joint_diskring_ksz_varweight_cmass_kendrick_pactf150daynight20200228maskgal60_cmass_kendrick_pactf90daynight20200228maskgal60_bootstrap.pdf ./cmass/

# summary tSZ
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/cmass_kendrick_tilecpactynocib/diskring_tsz_uniformweight.pdf ./cmass/
# cov tSZ
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/cmass_kendrick_tilecpactynocib/cor_diskring_tsz_uniformweight_bootstrap.pdf ./cmass/


# null and foreground tests
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/pipenulltests_ksz_150_cmass.pdf ./cmass/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/fgnulltests_ksz_150_cmass.pdf ./cmass/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/pipenulltests_tsz_150_cmass.pdf ./cmass/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/fgnulltests_tsz_150_cmass.pdf ./cmass/


# ksz bias from tsz, as a function of mMax
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/cmass_kendrick_pactf150daynight20200228maskgal60/diskring_ksz_varweight_mmax_tsztoksz.pdf ./cmass/


####################################################################################
# LOWZ

# summary ksz
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/lowz_kendrick_pactf150daynight20200228maskgal60/summary_ksz_150_90_cmass_kendrick.pdf ./lowz/
# joint cov 150-90
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/lowz_kendrick_pactf150daynight20200228maskgal60/cor_joint_diskring_ksz_varweight_cmass_kendrick_pactf150daynight20200228maskgal60_cmass_kendrick_pactf90daynight20200228maskgal60_bootstrap.pdf ./lowz/

# summary tSZ
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/lowz_kendrick_tilecpactynocib/diskring_tsz_uniformweight.pdf ./lowz/
# cov tSZ
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/lowz_kendrick_tilecpactynocib/cor_diskring_tsz_uniformweight_bootstrap.pdf ./lowz/


# null and foreground tests
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/pipenulltests_ksz_150_lowz.pdf ./lowz/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/fgnulltests_ksz_150_lowz.pdf ./lowz/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/pipenulltests_tsz_150_lowz.pdf ./lowz/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/fgnulltests_tsz_150_lowz.pdf ./lowz/


# ksz bias from tsz, as a function of mMax
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/lowz_kendrick_pactf150daynight20200228maskgal60/diskring_ksz_varweight_mmax_tsztoksz.pdf ./lowz/


####################################################################################
# pipeline tests

# signal mock maps
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/cmb_map/visu/* ./pipeline/

# pipeline end-to-end test and 2-halo terms
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/cmass_mariana_pactf150night20190311_test_endtoend_count_dirac_carmanu/test_mean_stacked_temperature_diskring_full.pdf ./pipeline/

# cov from GRF mocks 
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/compare_std_diskring_ksz_varweight_mocks0-800.pdf ./pipeline/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/compare_std_diskring_tsz_varweight_mocks0-800.pdf ./pipeline/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/cmass_mariana_mockgrf0_pactf150night20190311/cor_corTopMocksBottomBootstrap_diskring_ksz_varweight_mocks0-800.pdf ./pipeline/
rsync -avc cori:/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/cmass_mariana_mockgrf0_pactf150night20190311/cor_corTopMocksBottomBootstrap_diskring_tsz_varweight_mocks0-800.pdf ./pipeline/


