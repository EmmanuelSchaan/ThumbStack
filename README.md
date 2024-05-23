# ThumbStack: Moving lens branch

**Generating simulated maps of the moving lens effect:**

- flat sky: flatMLpainter/flatML.py
- curved sky: A corrected version of AstroPaint can be found here https://github.com/ali-beheshti/Astro-Paint

**Stacked estimators**

ThumbStack produces stacked maps and profiles, given catalogs of object positions and maps.
It is designed for thermal and kinematic Sunyaev-Zel'dovich measurements, and for moving lens measurements.
It outputs 2D stacked maps and radial profiles for different filters (e.g., aperture photometry filters), 
as well as their covariances, estimated through several methods including bootstrap.

**Citations**

ThumbStack was used in the following publication:
https://ui.adsabs.harvard.edu/abs/2020arXiv200905557S/abstract
Please feel free to use this code in your work.
If you use it in your publication, please cite the paper(s) above.
Contact Ali Beheshti and Emmanuel Schaan with any questions or suggestions (ali.beheshti@pitt.edu, eschaan@stanford.edu).

