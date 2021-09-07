# eMagLS
The End-to-End Magnitude Least Squares Binaural Renderer for Spherical Microphone Array Signals.

This repository contains MATLAB functions to obtain binaural rendering filters for the least squares (LS) method, the magnitude least squares (MagLS) method, and the newly proposed End-to-End MagLS methods using SH-domain processing (eMagLS) and using (raw) microphone signals (eMagLS2).
For more information and if you want to reference the code please refer to
   
   ```
   T. Deppisch, H. Helmholz, J. Ahrens, "End-to-End Magnitude Least Squares Binaural Rendering 
   of Spherical Microphone Array Signals," International 3D Audio Conference (I3DA), 2021.
   ```
   
The file `testEMagLS.m` contains an example on how to obtain the filters and apply them to an SMA recording. It also provides the opportunity to listen and compare the different renderers.

## Related Work
This repository depends on the [Spherical Harmonic Transform Library](https://github.com/polarch/Spherical-Harmonic-Transform) and uses a publicly accessible [HRIR set](https://zenodo.org/record/3928297) for demonstration purposes. It further uses an excerpt of an em32 recording from the [3D-MARCo library](https://zenodo.org/record/3477602). Files taken from external sources are subject to their corresponding licenses.

## Acknowledgment
We thank Facebook Reality Labs Research for funding this project.

## License
This software is licensed under a Non-Commercial Software License (see LICENSE for full details).
