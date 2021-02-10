#!/bin/bash


############################################################
# PACT

# generate full mask: footprint + point sources + Planck Milky Way
# same for 90 and 150GHz
#echo "generate masks PACT"
#python generate_mask_pact20200228_r2.py '60'
#python generate_mask_pact20200228_r2.py '70'

# reconvolve the 150GHz map to the beams of 
# the 90GHz, TileC, TileC deproj,
# for null tests
#python reconvolve_pact_r2.py

############################################################
# TileC

# generate the masks for TileC BOSS N and D56
#echo "generating masks"
#python generate_mask_tilec_v1.2.py 'cmbksz_d56'
#python generate_mask_tilec_v1.2.py 'cmbksz_boss'

# combine D56+BN into single maps
#echo "combining BN annd D56 maps"
#python combine_tilec_maps_v1.2.py "cmbksz"
#python combine_tilec_maps_v1.2.py "cmbksznoy"
#python combine_tilec_maps_v1.2.py "cmbksznocib"
#python combine_tilec_maps_v1.2.py "y"
#python combine_tilec_maps_v1.2.py "ynocib"
#python combine_tilec_maps_v1.2.py "ynocmb"


# reconvolve the TileC maps without deproj
# to the beam of the deproj maps,
# for null tests
#echo "Reconvolve TileC to TileC deproj, for null tests"
#python reconvolve_tilec_v1.2.py


############################################################
# Generate difference maps, for null tests
# this involves PACT maps and TileC maps

#python generate_diff_maps_r2.py


############################################################
# Make output readable by everyone

#chmod o+x ./output/thumbstack/*
#chmod o+r -R ./output/thumbstack/*

