# eMagLS
The End-to-End Magnitude Least Squares Binaural Renderer for Spherical Microphone Array Signals.

This repository contains MATLAB functions to obtain binaural rendering filters for the least squares (LS) method, the magnitude least squares (MagLS) method, and the newly proposed End-to-End MagLS methods using SH-domain processing (eMagLS) and using (raw) microphone signals (eMagLS2).
For more information and if you want to reference the code please refer to
   
   ```
   T. Deppisch, H. Helmholz, J. Ahrens, "End-to-End Magnitude Least Squares Binaural Rendering 
   of Spherical Microphone Array Signals," International Conference on Immersive and 3D Audio (I3DA), 2021.
   ```
   
The file `testEMagLS.m` contains an example on how to obtain the filters and apply them to an SMA recording. It also provides the opportunity to listen and compare the different renderers.

Make sure to clone the repository including submodules `git clone --recurse-submodules` or add the [Spherical Harmonic Transform Library](https://github.com/polarch/Spherical-Harmonic-Transform) manually to the `dependencies/` folder.

## Related Work
This repository depends on the [Spherical Harmonic Transform Library](https://github.com/polarch/Spherical-Harmonic-Transform) and uses a publicly accessible [HRIR set](https://zenodo.org/record/3928297) for demonstration purposes. It further uses an excerpt of an em32 recording from the [3D-MARCo library](https://zenodo.org/record/3477602). Files from external sources are subject to their corresponding licenses.

## Changelog
### 2022-05-18
- Update `testEMagLs.m` to use a different SH basis implementation more easily
- Update `getRadialFilter.m` and `getSMAIRMatrix.m` to use `sphModalCoeffs()` from Array-Response-Simulator</br>
(therefore also remove own implementations of spherical hankel and bessel functions)</br>
(this causes the resulting eMagLS and eMagLS2 rendering filters to be slightly different at 0 Hz and the Nyquist bin)</br>
(therefore the reference for eMagLS and eMagLS2 rendering filters are updated for verification)
- Update `getSMAIRMatrix.m` to streamline conjugate and transpose operations
- Update `getEMagLs2Filters.m` to streamline conjugate and transpose operations (results in neglectable maximum normalized absolute difference of 1e-14)
- Update `getEMagLsFilters.m` to streamline conjugate and transpose operations (results in neglectable maximum normalized absolute difference of 1e-15)
- Update `getLsFilters.m` to streamline conjugate and transpose operations (results in neglectable maximum normalized absolute difference of 1e-16)
- Update `getMagLsFilters.m` to streamline conjugate and transpose operations (results in neglectable maximum normalized absolute difference of 1e-15)
- Update `testEMagLs.m` to accept tolerance when verifying rendering filters against the provided reference
- Update functions to use `ifft()` without the forced symmetric parameter (fix to yield real or complex filters depending on the basis type)
- Update `getEMagLsFilters.m` and `getEMagLs2Filters.m` to slightly improve computation time
### 2022-05-03
- Update `testEMagLs.m` to add used SH basis type in exported file names (update file names of respective reference results)
- Update `testEMagLs.m` to export audio rendering results as an option
- Fix functions to use case-insensitive string comparisons for parameters
- Update functions to provide SH basis implementation function as an optional parameter (see `testEMagLs.m` for a usage example)
- Update functions to provide utilized constants at the top
- Update functions to provide SH basis type (according to basis implementation) as an optional parameter
- Update `testEMagLs.m` to add verification of rendering filters against the provided reference
- Update `testEMagLs.m` to add global configuration variables
- Update `testEMagLs.m` to export reference results for LS, MagLS, eMagLS and eMagLS2 filters (also add the respective reference results)
- Update `testEMagLs.m` to provide configuration of the MagLS diffuseness constraint
- Update `testEMagLs.m` to be more verbose (reintroduce length limitation of the rendered SMA recording)
- Remove default SMA recording dataset (update `testEMagLs.m` to download data on demand)
- Remove default HRIR dataset (update `testEMagLs.m`t to download data on demand)
- Fix all function endings and formatting
### 2022-03-16
- Update `testEMagLs.m` to be verbose about audio playback of rendering results
- Update `README.md` with information on initializing git submodules 
### 2021-10-10
- Initial code release

## Acknowledgment
We thank Facebook Reality Labs Research for funding this project.

## License
This software is licensed under a Non-Commercial Software License (see [LICENSE](https://github.com/thomasdeppisch/eMagLS/blob/master/LICENSE) for full details).
