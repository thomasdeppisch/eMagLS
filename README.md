# End-to-End Magnitude Least Squares Binaural Rendering
This repository contains MATLAB functions to obtain binaural rendering filters for the least squares (LS) method, the magnitude least squares (MagLS) method, and the newly proposed End-to-End MagLS methods using SH-domain processing (eMagLS) and using (raw) microphone signals (eMagLS2).

The repository has been extended to provide rendering filters for the signals from spherical microphone arrays (SMAs) and equatorial microphone arrays (EMAs).

For more information and if you want to reference the code please refer to the [publication for SMAs](https://research.chalmers.se/publication/528436/file/528436_Fulltext.pdf) or the [publication for EMAs](https://research.chalmers.se/publication/535525/file/535525_Fulltext.pdf) respectively.
   
   ```
   T. Deppisch, H. Helmholz, and J. Ahrens,
   “End-to-End Magnitude Least Squares Binaural Rendering of Spherical Microphone Array Signals,”
   in 2021 Immersive and 3D Audio: from Architecture to Automotive (I3DA), 2021, pp. 1–7. doi: 10.1109/I3DA48870.2021.9610864.
   ```
   ```
   H. Helmholz, T. Deppisch, and J. Ahrens,
   “End-to-End Magnitude Least Squares Binaural Rendering for Equatorial Microphone Arrays,”
   in Fortschritte der Akustik -- DAGA 2023, 2023, pp. 1679–1682.
   ```
   
The file [testEMagLs.m](testEMagLs.m) contains an example on how to obtain the filters and apply them to an SMA recording. It also provides the opportunity to listen and compare the different renderers.

Make sure to clone the repository including submodules `git clone --recurse-submodules` or add the [Spherical Harmonic Transform Library](https://github.com/polarch/Spherical-Harmonic-Transform) and the [Ambisonic Encoding Toolbox](https://github.com/AppliedAcousticsChalmers/ambisonic-encoding) manually to the `dependencies/` folder.

## Related Work
This repository depends on the [Spherical Harmonic Transform Library](https://github.com/polarch/Spherical-Harmonic-Transform), the [Array Response Simulator](https://github.com/polarch/Array-Response-Simulator.git), the [Ambisonic Encoding Toolbox](https://github.com/AppliedAcousticsChalmers/ambisonic-encoding) and uses a publicly accessible [HRIR set](https://zenodo.org/record/3928297) for demonstration purposes. It further uses an excerpt of an em32 recording from the [3D-MARCo library](https://zenodo.org/record/3477602). Files from external sources are subject to their corresponding licenses.

## Acknowledgment
We thank Facebook Reality Labs Research for funding this project.

## Changelog
See [CHANGELOG](CHANGELOG.md) for full details.

## License
This software is licensed under a Non-Commercial Software License (see [LICENSE](LICENSE) for full details).
