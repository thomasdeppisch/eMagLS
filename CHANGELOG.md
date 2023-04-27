# Changelog

### Unreleased
- Update `getRadialFilter.m` to set the Nyquist bin for even-length filters after calculation
- Update eMagLS functions with error message when applying the diffuseness constraint with "complex" SH conventions (this is not implemented yet)
- Update `testEMagLs.m` to not perform validation of the rendering filters against the reference by default
- Update eMagLS functions documentation when manually setting the DC bin after calculation
- Update eMagLS functions to manually set the DC bin after calculation</br>
(therefore the reference for the rendering filters are updated for verification)
- Update `testEMagLs.m` to download the MIRO class if it is not available for loading the default HRIR data set
- Update `testEMagLs.m` documentation
- Fix missing `Spherical-Harmonic-Transform` git submodule
- Update header documentation in filter functions
- Add `CHANGELOG.md` with changelog moved from `README.md`
- Update `README.md` with references to the EMA publication 
- Add `ambisonic-encoding` as git submodule (required for `get_sound_field_sh_coeffs_from_ema_t()` in `getEMagLsFiltersEMAinSH.m`)
- Update header documentation of all functions
- Remove outdated functions that were used only in preliminary implementations of EMA rendering filters
- Update `applyRadialFilter.m` to perform causalization and fading analogous to other filtering functions
- Resolve merge conflicts to main branch
### 2023-04-20
- Update `getMagLsFilters.m` to improve computation speed by matricizing all operations for both ears</br>
(something similar should also be implemented for the other functions in the future)</br>
(this causes the resulting rendering filters with diffuseness constraint to be different at high frequencies in a non-meaningful way)</br>
(therefore the reference for the rendering filters with diffuseness constraint are updated for verification)
- Remove outdated functions with preliminary implementations of EMA rendering filters
- Update eMagLS functions with comment that diffuseness constraint causes slight magnitude deviations at low frequencies
- Fix EMA functions to implement the diffuseness constraint
- Update SMA functions to improve syntax for diffuseness constraint (the resulting filters are identical)
- Add reference rendering filters with diffuseness constraint for verification
- Update `testEMagLs.m` to yield reference names including the status of the diffuseness constraint</br>
(rename the respective reference rendering filter files)
- Fix `getEMagLs2Filters.m` when applying the diffuseness constraint
### 2023-03-31
- Update functions to use a transition frequency of at least 1 kHz (this makes a difference only for SH order 1)
- Update `getEMagLsFiltersEMAinSH.m` to improve computation speed by performing all sound field rotations in one operation
- Fix `getEMagLsFiltersEMAinSH.m` vertical rotation direction (it was inverted before)
- Update `getSMAIRMatrix.m` to work for higher requested array orders than estimated simulation order
- Update functions to improve computation speed by pre-allocating complex matrices when needed
- Update `getEMagLsFiltersEMAinSH.m` to use existing EMA function from the Ambisonic Encoding toolbox
- Update `getMagLsSphericalHeadFilter.m` and `getMagLsArrayDiffuseFilter.m` to remove amplitude normalization and yield no magnitude changes at low frequencies
### 2023-03-08
- Add `getEMagLsFiltersEMAinSH.m` to provide the implementation of EMA rendering filters in spherical harmonics
- Update functions to calculate the HRIR group delay in a simpler way without meaningful changes in the result</br>
(this causes the resulting MagLS, eMagLS and eMagLS2 rendering filters to be different in a non-meaningful way at very high frequencies)</br>
(therefore the reference for MagLS, eMagLS and eMagLS2 rendering filters are updated for verification)
- Add functions `getMagLsSphericalHeadFilter.m` and `getMagLsArrayDiffuseFilter.m` to provide direction-independent array equalization filters that can be useful in combination with MagLS
- Update `getEMagLsFiltersEMAinCH.m`, `getEMagLsFiltersEMAinCHtoSH.m` and `getEMagLsFiltersEMAinEHtoSH.m` to improve computation speed</br>
(extract feasible matrix multiplications to be pre-computed before the loop)
- Update `getEMagLsFiltersEMAinCHtoSH.m` and `getEMagLsFiltersEMAinEHtoSH.m` to extract helper functions</br>
(add the respective `ch_fromShIds.m` and `eh_fromShIds.m` functions)
- Fix `getEMagLsFiltersEMAinEHtoSH.m` to reflect switch from "soundfieldsynthesis" scripts to updated "Ambisonic Encoding" toolbox
- Fix functions to calculate the filters at the specified oversampled length</br>
(now the oversampling actually uses double the specified length up to a maximum of 2048 samples)</br>
(this causes all the resulting rendering filters to be different at high frequencies in a non-meaningful way)</br>
(therefore the reference for all the rendering filters are updated for verification)
### 2022-11-18
- Add `getEMagLsFiltersEMAinCH.m` and `getEMagLsFiltersEMAinCHtoSH.m` with different preliminary implementations of EMA rendering filters
- Rename `getEMagLsFiltersEMA.m` into ``getEMagLsFiltersEMAinEHtoSH.m` to reflect its function
- Fix `getEMagLsFiltersEMA.m` parsing of default CH basis function handle
- Update `getEMagLsFiltersEMA.m` SVD regularization parameter to yield more sensible rendering filters that</br>
(better behaviour in time domain at the cost of small deviations in magnitude response)
- Update `getEMagLsFiltersEMA.m` to temporarily use only horizontal HRIRs for magnitude optimization</br>
(spherical HRIR data will be SH sub-sampled to a horizontal grid)
- Update functions to simplify syntax of singular value decomposition
- Update functions to extract generation of time domain window for resulting rendering filters</br>
(add the respective `getFadeWindow.m` function)
### 2022-09-19
- Update `getEMagLsFiltersEMA.m` to remove unused array grid zenith coordinates</br>
(add non-functioning example call to `testEMagLs.m`)
- Add `getEMagLsFiltersEMA.m` to generate eMagLS rendering filters for equatorial microphone arrays (EMA)</br>
(add the respective `getCH.m` and `ch_stackOrder.m` functions)</br>
(this is does not have an example configuration in `testEMagLs.m` yet)</br>
(the rendering filters for "complex" SH basis types do not function as intended yet)</br>
(the implementation of the diffuseness constraint has not been verified yet)
- Add reference rendering filters for "complex" SH convention for verification
- Update `getMagLsFilters.m`, `getEMagLsFilters.m` and `getEMagLs2Filters.m` to use delays with subsample precision and restore original inter-aural group delay difference after magnitude least-squares optimization</br>
(add the respective `applySubsampleDelay.m` function)</br>
(this causes the resulting MagLS, eMagLS and eMagLS2 rendering filters to be different in a non-meaningful way at very high frequencies)</br>
(therefore the reference for MagLS, eMagLS and eMagLS2 rendering filters are updated for verification)
### 2022-06-15
- Update `getSMAIRMatrix.m` to determine the simulation order based on the configuration aliasing frequency</br>
(this causes the resulting eMagLS and eMagLS2 rendering filters to be different in a non-meaningful way at very high frequencies)</br>
(therefore the reference for eMagLS and eMagLS2 rendering filters are updated for verification)
- Update `getMagLsFilters.m`, `getEMagLsFilters.m` and `getEMagLs2Filters.m` to use a transition frequency of 500 Hz times SH order instead of always 2000 Hz</br>
(changes the signature of `getEMagLs2Filters.m` to require the desired SH order)</br>
(yields no changes for the example Eigenmike configuration)
- Update `testEMagLs.m` to export binaural renderings as WAV with 64bit resolution
- Update functions to work for "complex" SH conventions by computing double-sided spectra</br>
(required in `getMagLsFilters.m` and `getEMagLsFilters.m` but not in `getEMagLs2Filters.m`)
- Update `testEMagLs.m` to alternatively evaluate spectral difference in verification of rendering filters against provided reference
- Update functions to compute 0 Hz bin not separately</br>
(this causes the resulting eMagLS and eMagLS2 rendering filters to be different)</br>
(therefore the reference for eMagLS and eMagLS2 rendering filters are updated for verification)
- Update `testEMagLs.m` to verify the SH convention "wikipedia" against the "real" reference
### 2022-05-19
- Update `binauralDecode.m` to warn when rendering discards imaginary signal parts (may occur for "complex" SH basis functions)
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
- Update functions to use `ifft()` without the forced symmetric parameter (fix to yield "real" or "complex" filters depending on the basis type)
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
