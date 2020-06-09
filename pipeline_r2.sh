#!/bin/bash


############################################################
# PACT

# generate full mask: footprint + point sources + Planck Milky Way
# same for 90 and 150GHz
#python generate_mask_pact20200228_r2.py '60'
#python generate_mask_pact20200228_r2.py '70'

# reconvolve the 150GHz map to the beam of the 90GHz, 
# for null tests
#python reconvolve_pact_150to90_r2.py

############################################################
# TileC

# generate the masks for TileC BOSS N and D56
#python generate_mask_tilec_v1.2.py 'cmbksz_d56'
#python generate_mask_tilec_v1.2.py 'cmbksz_boss'

# combine D56+BN into single maps
#python combine_tilec_maps_v1.2.py "cmbksz"
#python combine_tilec_maps_v1.2.py "cmbksznoy"
#python combine_tilec_maps_v1.2.py "cmbksznocib"
#python combine_tilec_maps_v1.2.py "y"
#python combine_tilec_maps_v1.2.py "ynocib"
#python combine_tilec_maps_v1.2.py "ynocmb"

# Reconvolve TileC maps to 1.4'
#python reconvolve_tilec_v1.2.py "cmbksz"
#python reconvolve_tilec_v1.2.py "cmbksz_d56"
#python reconvolve_tilec_v1.2.py "cmbksz_boss"
#python reconvolve_tilec_v1.2.py "cmbksznoy"
#python reconvolve_tilec_v1.2.py "cmbksznoy_d56"
#python reconvolve_tilec_v1.2.py "cmbksznoy_boss"
#python reconvolve_tilec_v1.2.py "cmbksznocib"
#python reconvolve_tilec_v1.2.py "cmbksznocib_d56"
#python reconvolve_tilec_v1.2.py "cmbksznoy_boss"
#python reconvolve_tilec_v1.2.py "y"
#python reconvolve_tilec_v1.2.py "y_d56"
#python reconvolve_tilec_v1.2.py "y_boss"
#python reconvolve_tilec_v1.2.py "ynocib"
#python reconvolve_tilec_v1.2.py "ynocib_d56"
#python reconvolve_tilec_v1.2.py "ynocib_boss"
#python reconvolve_tilec_v1.2.py "ynocmb"
#python reconvolve_tilec_v1.2.py "ynocmb_d56"
#python reconvolve_tilec_v1.2.py "ynocmb_boss"


############################################################
# Generate difference maps, for null tests
# this involves PACT maps and TileC maps

#python generate_diff_maps_r2.py


############################################################
# Make output readable by everyone

#chmod o+x ./output/thumbstack/*
#chmod o+r -R ./output/thumbstack/*

