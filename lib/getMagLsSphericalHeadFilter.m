function [wShf, W_Shf] = getMagLsSphericalHeadFilter(micRadius, order, fs, len)
% [wShf, W_Shf] = getMagLsSphericalHeadFilter(micRadius, order, fs, len)
%
% This function calculates the inverted spherical head filter to be 
% optionally applied in addition to the MagLS rendering filters.
% For more information, please refer to
%   Ben-Hur, Brinkmann, Sheaffer, Weinzierl, Rafaely,
%   "Spectral Equalization in Binaural Signals Represented by Order-Truncated Spherical Harmonics",
%   J. Acoust. Soc. Am., vol. 141, no. 6, pp. 4087–4096, 2017, doi: 10.1121/1.4983652.
% 
% wShf                   .. time-domain filter (linear-phase-like with fade-out window)
% W_Shf                  .. frequency-domain filter (zero-phase-like)
% micRadius              .. array radius in m
% order                  .. array SH order
% fs                     .. sampling frequency in Hz
% len                    .. desired length of filter
% 
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Hannes Helmholz, 2023

NFFT_MAX_LEN            = 2048; % maxium oversamping length in samples
SIMULATION_ARRAY_TYPE   = 'rigid'; % see `sphModalCoeffs()`
SIMULATION_MIC_DIR      = 0; % see `sphModalCoeffs()`
SIMULATION_C            = 343; % speed of sound in m/s

nfft = min(NFFT_MAX_LEN, 2 * len); % apply frequency-domain oversampling
f = linspace(0, fs/2, nfft/2+1).';
kr = 2*pi*f/SIMULATION_C * micRadius;
simulationOrder = ceil(fs * pi * micRadius / SIMULATION_C);

% get modal radial filters
bn_Hi = sphModalCoeffs(simulationOrder, kr, SIMULATION_ARRAY_TYPE, SIMULATION_MIC_DIR);
bn_Lo = bn_Hi(:, 1:order+1);

% expand to SHs
bn_Hi = sh_repToOrder(bn_Hi.').';
bn_Lo = sh_repToOrder(bn_Lo.').';

% calculate diffuse-field response
bn_Hi_df = rms(abs(bn_Hi), 2) * sqrt(size(bn_Hi, 2)) / (4*pi);
bn_Lo_df = rms(abs(bn_Lo), 2) * sqrt(size(bn_Lo, 2)) / (4*pi);

% calculate scattering diffuse-field difference (Spherical Head Filter)
W_Shf = bn_Hi_df ./ bn_Lo_df;

% invert spectrum
W_Shf = 1 ./ W_Shf;

% transform into time domain
numPosFreqs = nfft/2+1;
W_Shf = [W_Shf(1:numPosFreqs, :); flipud(conj(W_Shf(2:numPosFreqs-1, :)))];
wShf = ifft(W_Shf);

% shift from zero-phase-like to linear-phase-like
n_shift = nfft/2;
wShf = applySubsampleDelay(wShf, n_shift);

% shorten to target length
wShf = wShf(n_shift-len/2+1:n_shift+len/2, :);

% fade
fade_win = getFadeWindow(len);
wShf = wShf .* fade_win;

end
