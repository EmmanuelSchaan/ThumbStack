#!/bin/bash


# combine D56+BN into single maps
#python combine_tilec_maps.py "cmbksz"
#python combine_tilec_maps.py "cmbksznoy"
#python combine_tilec_maps.py "cmbksznocib"
#python combine_tilec_maps.py "y"
#python combine_tilec_maps.py "ynocib"



## Reconvolve TileC maps to 1.4'
python reconvolve_tilec.py "cmbksz"
python reconvolve_tilec.py "cmbksz_d56"
python reconvolve_tilec.py "cmbksz_boss"
python reconvolve_tilec.py "cmbksznoy"
python reconvolve_tilec.py "cmbksznoy_d56"
python reconvolve_tilec.py "cmbksznoy_boss"
#python reconvolve_tilec.py "cmbksznocib"
#python reconvolve_tilec.py "cmbksznocib_d56"
#python reconvolve_tilec.py "cmbksznoy_boss"
python reconvolve_tilec.py "y"
python reconvolve_tilec.py "y_d56"
python reconvolve_tilec.py "y_boss"
#python reconvolve_tilec.py "ynocib"
#python reconvolve_tilec.py "ynocib_d56"
#python reconvolve_tilec.py "ynocib_boss"

