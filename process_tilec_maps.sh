#!/bin/bash



# generate the masks for TileC BOSS N and D56
#python read_map_tilec_pact_v1.2.py 'cmbksz_d56'
#python read_map_tilec_pact_v1.2.py 'cmbksz_boss'


# combine D56+BN into single maps
#python combine_tilec_maps_v1.2.py "cmbksz"
#python combine_tilec_maps_v1.2.py "cmbksznoy"
#python combine_tilec_maps_v1.2.py "cmbksznocib"
#python combine_tilec_maps_v1.2.py "y"
#python combine_tilec_maps_v1.2.py "ynocib"



# Reconvolve TileC maps to 1.4'
python reconvolve_tilec_v1.2.py "cmbksz"
python reconvolve_tilec_v1.2.py "cmbksz_d56"
python reconvolve_tilec_v1.2.py "cmbksz_boss"
python reconvolve_tilec_v1.2.py "cmbksznoy"
python reconvolve_tilec_v1.2.py "cmbksznoy_d56"
python reconvolve_tilec_v1.2.py "cmbksznoy_boss"
python reconvolve_tilec_v1.2.py "cmbksznocib"
python reconvolve_tilec_v1.2.py "cmbksznocib_d56"
python reconvolve_tilec_v1.2.py "cmbksznoy_boss"
python reconvolve_tilec_v1.2.py "y"
python reconvolve_tilec_v1.2.py "y_d56"
python reconvolve_tilec_v1.2.py "y_boss"
python reconvolve_tilec_v1.2.py "ynocib"
python reconvolve_tilec_v1.2.py "ynocib_d56"
python reconvolve_tilec_v1.2.py "ynocib_boss"

