function [wMlsL, wMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, micGridZenRad, order, fs, len, shDefinition, shFunction)
% [wMlsL, wMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
%     micRadius, micGridAziRad, micGridZenRad, order, fs, len, shDefinition, shFunction)
%
% This function calculates eMagLS binaural decoding filters for spherical microphone arrays.
% For more information, please refer to
%   T. Deppisch, H. Helmholz, and J. Ahrens,
%   “End-to-End Magnitude Least Squares Binaural Rendering of Spherical Microphone Array Signals,”
%   in 2021 Immersive and 3D Audio: from Architecture to Automotive (I3DA), 2021, pp. 1–7. doi: 10.1109/I3DA48870.2021.9610864.
%
% wMlsL                  .. time-domain decoding filter for left ear
% wMlsR                  .. time-domain decoding filter for right ear
% hL                     .. HRIR set for left ear (numSamples x numDirections)
% hR                     .. HRIR set for right ear (numSamples x numDirections)
% hrirGridAziRad         .. grid azimuth angles in radians of HRIR set (numDirections x 1)
% hrirGridZenRad         .. grid zenith angles in radians of HRIR set (numDirections x 1)
% micRadius              .. radius of SMA
% micGridAziRad          .. SMA grid azimuth angles in radians
% micGridZenRad          .. SMA grid zenith angles in radians
% order                  .. SH output order
% fs                     .. sampling frequency in Hz
% len                    .. desired length of eMagLS filters
% shDefinition           .. SH basis type according to utilized shFunction, default: 'real'
% shFunction             .. SH basis function (see testEMagLs.m for example), default: @getSH
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Thomas Deppisch, 2023

if nargin < 12; shFunction = @getSH; end
if nargin < 11 || isempty(shDefinition); shDefinition = 'real'; end

NFFT_MAX_LEN            = 2048; % maxium oversamping length in samples
F_CUT_MIN_FREQ          = 1e3; % minimum transition freqeuncy in Hz
SIMULATION_WAVE_MODEL   = 'planeWave'; % see `getSMAIRMatrix()`
SIMULATION_ARRAY_TYPE   = 'rigid'; % see `getSMAIRMatrix()`
SVD_REGUL_CONST         = 0.01;

% TODO: Implement dealing with HRIRs that are longer than the requested filter
assert(len >= size(hL, 1), 'len too short');

nfft = min(NFFT_MAX_LEN, 2 * len); % apply frequency-domain oversampling
f = linspace(0, fs/2, nfft/2+1).';
numPosFreqs = length(f);
f_cut = max(F_CUT_MIN_FREQ, 500 * order); % from N > k
k_cut = ceil(f_cut / f(2));

% simulate plane wave impinging on SMA
params.order = order;
params.fs = fs;
params.irLen = nfft;
params.oversamplingFactor = 1;
params.simulateAliasing = true;
params.radialFilter = 'none';
params.smaRadius = micRadius;
params.smaDesignAziZenRad = [micGridAziRad, micGridZenRad];
params.waveModel = SIMULATION_WAVE_MODEL;
params.arrayType = SIMULATION_ARRAY_TYPE;
params.shDefinition = shDefinition;
params.shFunction = shFunction;
smairMat = getSMAIRMatrix(params);
simulationOrder = sqrt(size(smairMat, 2)) - 1;

numHarmonics = (order+1)^2;
numDirections = size(hL, 2);
Y_Hi_conj = shFunction(simulationOrder, [hrirGridAziRad, hrirGridZenRad], shDefinition)';

% zero pad and remove group delay with subsample precision
% (alternative to applying global phase delay later)
hL(end+1:nfft, :) = 0;
hR(end+1:nfft, :) = 0;
grpDL = median(grpdelay(sum(hL, 2), 1, f, fs));
grpDR = median(grpdelay(sum(hR, 2), 1, f, fs));
hL = applySubsampleDelay(hL, -grpDL);
hR = applySubsampleDelay(hR, -grpDR);

% transform into frequency domain
HL = fft(hL);
HR = fft(hR);

W_MLS_l = zeros(nfft, numHarmonics, 'like', HL);
W_MLS_r = zeros(nfft, numHarmonics, 'like', HL);
for k = 2:numPosFreqs
    % positive frequencies
    pwGrid = smairMat(:,:,k) * Y_Hi_conj;
    [U, s, V] = svd(pwGrid.', 'econ', 'vector');
    s = 1 ./ max(s, SVD_REGUL_CONST * max(s)); % regularize
    Y_reg_inv = conj(U) * (s .* V.');

    if k < k_cut % least-squares below cut
        W_MLS_l(k,:) = HL(k,:) * Y_reg_inv;
        W_MLS_r(k,:) = HR(k,:) * Y_reg_inv;
    else % magnitude least-squares above cut
        phi_l = angle(W_MLS_l(k-1,:) * pwGrid);
        phi_r = angle(W_MLS_r(k-1,:) * pwGrid);
        if k == numPosFreqs % Nyquist bin, is even
            W_MLS_l(k,:) = real(abs(HL(k,:)) .* exp(1i * phi_l)) * Y_reg_inv;
            W_MLS_r(k,:) = real(abs(HR(k,:)) .* exp(1i * phi_r)) * Y_reg_inv;
        else
            W_MLS_l(k,:) = abs(HL(k,:)) .* exp(1i * phi_l) * Y_reg_inv;
            W_MLS_r(k,:) = abs(HR(k,:)) .* exp(1i * phi_r) * Y_reg_inv;
        end
    end
end

% mamnually set the DC bin (use `real()` instead of `abs()`, which causes
% strong magnitude errors in the rendering results at low frequencies)
W_MLS_l(1, :) = real(W_MLS_l(2, :));
W_MLS_r(1, :) = real(W_MLS_r(2, :));

% to time domain
if isreal(Y_Hi_conj)
    W_MLS_l = [W_MLS_l(1:numPosFreqs, :); flipud(conj(W_MLS_l(2:numPosFreqs-1, :)))];
    W_MLS_r = [W_MLS_r(1:numPosFreqs, :); flipud(conj(W_MLS_r(2:numPosFreqs-1, :)))];
else
    W_MLS_l = getShFreqDomainConjugate(W_MLS_l(1:numPosFreqs, :));
    W_MLS_r = getShFreqDomainConjugate(W_MLS_r(1:numPosFreqs, :));
end

wMlsL = ifft(W_MLS_l);
wMlsR = ifft(W_MLS_r);
if isreal(Y_Hi_conj)
    assert(isreal(wMlsL), 'Resulting decoding filters are not real valued.');
    assert(isreal(wMlsR), 'Resulting decoding filters are not real valued.');
end

% shift from zero-phase-like to linear-phase-like
% and restore initial group-delay difference between ears
n_shift = nfft/2;
wMlsL = applySubsampleDelay(wMlsL, n_shift);
wMlsR = applySubsampleDelay(wMlsR, n_shift+grpDR-grpDL);

% shorten to target length
wMlsL = wMlsL(n_shift-len/2+1:n_shift+len/2, :);
wMlsR = wMlsR(n_shift-len/2+1:n_shift+len/2, :);

% fade
fade_win = getFadeWindow(len);
wMlsL = wMlsL .* fade_win;
wMlsR = wMlsR .* fade_win;

end
