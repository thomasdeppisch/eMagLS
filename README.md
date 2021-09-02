# eMagLS
The End-to-End Magnitude Least Squares Binaural Renderer for Spherical Microphone Array Signals.

This repository contains MATLAB functions to obtain binaural rendering filters for the least squares (LS) method, the magnitude least squares (MagLS) method, and the newly proposed End-to-End MagLS methods using SH-domain processing (eMagLS) and using (raw) microphone signals (eMagLS2).
For more information and if you want to reference the code please refer to
   
   ```
   T. Deppisch, H. Helmholz, J. Ahrens, "End-to-End Magnitude Least Squares Binaural Rendering 
   of Spherical Microphone Array Signals," International 3D Audio Conference (I3DA), 2021.
   ```
   
The file `testEMagLS.m` contains an example on how to obtain the filters and apply them to a SMA recording. It also provides the opportunity to listen and compare the different renders.
