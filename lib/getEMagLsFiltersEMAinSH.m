function [wMlsL, wMlsR] = getEMagLsFiltersEMAinSH(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, order, fs, len, shDefinition, shFunction, chFunction)
% [wMlsL, wMlsR] = getEMagLsFiltersEMAinSH(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
%     micRadius, micGridAziRad, order, fs, len, shDefinition, shFunction)
%
% This function calculates eMagLS binaural decoding filters in spherical harmonics
% for equatorial microphone arrays, allowing for 3-DOF head rotations.
% For more information, please refer to
%   H. Helmholz, T. Deppisch, and J. Ahrens,
%   “End-to-End Magnitude Least Squares Binaural Rendering for Equatorial Microphone Arrays,”
%   in Fortschritte der Akustik -- DAGA 2023, 2023, pp. 1679–1682.
%
% wMlsL                  .. time-domain decoding filter for left ear
% wMlsR                  .. time-domain decoding filter for right ear
% hL                     .. HRIR set for left ear (numSamples x numDirections)
% hR                     .. HRIR set for right ear (numSamples x numDirections)
% hrirGridAziRad         .. grid azimuth angles in radians of HRIR set (numDirections x 1)
% hrirGridZenRad         .. grid zenith angles in radians of HRIR set (numDirections x 1)
% micRadius              .. radius of EMA
% micGridAziRad          .. EMA grid azimuth angles in radians
% order                  .. SH output order
% fs                     .. sampling frequency in Hz
% len                    .. desired length of eMagLS filters
% shDefinition           .. SH basis type according to utilized shFunction, default: 'real'
% shFunction             .. SH basis function (see testEMagLs.m for example), default: @getSH
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Hannes Helmholz, 2023

if nargin < 12; chFunction = @getCH; end
if nargin < 11; shFunction = @getSH; end
if nargin < 10 || isempty(shDefinition); shDefinition = 'real'; end

NFFT_MAX_LEN            = 2048; % maxium oversamping length in samples
F_CUT_MIN_FREQ          = 1e3; % minimum transition freqeuncy in Hz
SIMULATION_WAVE_MODEL   = 'planeWave'; % see `getSMAIRMatrix()`
SIMULATION_ARRAY_TYPE   = 'rigid'; % see `getSMAIRMatrix()`
SVD_REGUL_CONST         = 0.01;

% TODO: Implement dealing with HRIRs that are longer than the requested filter
assert(len >= size(hL, 1), 'len too short');

nfft = min(NFFT_MAX_LEN, 2 * len); % apply frequency-domain oversampling
f = linspace(0, fs/2, nfft/2+1).';
f_cut = max(F_CUT_MIN_FREQ, 500 * order); % from N > k
k_cut = ceil(f_cut / f(2));

% simulate plane wave impinging on EMA
params.returnRawMicSigs = true; % raw mic signals, no SHs!
params.order = order;
params.fs = fs;
params.irLen = nfft;
params.oversamplingFactor = 1;
params.radialFilter = 'none';
params.smaRadius = micRadius;
params.smaDesignAziZenRad = [micGridAziRad, ones(size(micGridAziRad)) * pi/2];
params.waveModel = SIMULATION_WAVE_MODEL;
params.arrayType = SIMULATION_ARRAY_TYPE;
params.shDefinition = shDefinition;
params.shFunction = shFunction;
emaIrMat = getSMAIRMatrix(params);
simulationOrder = sqrt(size(emaIrMat, 2)) - 1;

% sample the microphone responses to PW directions on the EMA from all HRIR
% directions, but mapped to elevation = 0
Y_hor_conj = shFunction(simulationOrder, [hrirGridAziRad, ones(size(hrirGridAziRad)) * pi/2], shDefinition)';
emaIrDir = pagemtimes(emaIrMat, Y_hor_conj);
clear emaIrMat;

numPosFreqs = length(f);
emaIrDir = permute(emaIrDir, [3, 1, 2]);

% decompose with EMA functions without radial filters
numHarmonics = (order+1)^2;
numDirections = size(hL, 2);
emaIrDir_sh = zeros(numPosFreqs, numHarmonics, numDirections);
YCh = chFunction(order, micGridAziRad, shDefinition);
J = getChToShExpansionMatrix(order,shDefinition);
for d = 1 : numDirections
    emaIrDir_sh(:, :, d) = emaIrDir(:, :, d) * pinv(YCh.') * J.';
end

% rotate the EMA SH sound field to impose the HRIR elevation
for d = 1 : numDirections
    if hrirGridZenRad(d) ~= pi/2
        % NOTE: This uses the SHT convention of spherical harmonics.
        %       Therefore, the resulting rotation matrices may not be
        %       correct for other conventions!

        % 1) rotate SH represention so that the plane wave impinges from the front
        % 2) rotate SH represention to impose the inverted desired elevation of incidence
        %    (and convert colatitude to elevation)
        % 3) rotate SH represention back to the original azimuth
        euler_matrix  = euler2rotationMatrix(-hrirGridAziRad(d), ...
            hrirGridZenRad(d) - pi/2, hrirGridAziRad(d), 'zyz');
        sh_rot_matrix = getSHrotMtx(euler_matrix, order, shDefinition);
        emaIrDir_sh(:, :, d) = emaIrDir_sh(:, :, d) * sh_rot_matrix;
    end
end; clear d euler_matrix sh_rot_matrix;

% to frequency domain
pwGridAll = permute(emaIrDir_sh, [2, 3, 1]);
clear emaIrDir_sh;

% zero pad and remove group delay with subsample precision
% (alternative to applying global phase delay later)
hL(end+1:nfft, :) = 0;
hR(end+1:nfft, :) = 0;
grpDL = median(grpdelay(sum(hL, 2), 1, f, fs));
grpDR = median(grpdelay(sum(hR, 2), 1, f, fs));
hL = applySubsampleDelay(hL, -grpDL);
hR = applySubsampleDelay(hR, -grpDR);

% to frequency domain
HL = fft(hL,nfft);
HR = fft(hR,nfft);

W_MLS_l = zeros(nfft, numHarmonics, 'like', HL);
W_MLS_r = zeros(nfft, numHarmonics, 'like', HL);
for k = 2:numPosFreqs
    % positive frequencies
    pwGrid = pwGridAll(:,:,k);
    [U, s, V] = svd(pwGrid.', 'econ', 'vector');
    s = 1 ./ max(s, SVD_REGUL_CONST * max(s)); % regularize
    Y_reg_inv = conj(U) * (s .* V.');

    if k < k_cut % least-squares below cut
        W_MLS_l(k,:) = HL(k,:) * Y_reg_inv;
        W_MLS_r(k,:) = HR(k,:) * Y_reg_inv;
    else % magnitude least-squares above cut
        phi_l = angle(W_MLS_l(k-1,:) * pwGrid);
        phi_r = angle(W_MLS_r(k-1,:) * pwGrid);
        if k == numPosFreqs % Nyquist bin
            W_MLS_l(k,:) = real(abs(HL(k,:)) .* exp(1i * phi_l)) * Y_reg_inv;
            W_MLS_r(k,:) = real(abs(HR(k,:)) .* exp(1i * phi_r)) * Y_reg_inv;
        else
            W_MLS_l(k,:) = abs(HL(k,:)) .* exp(1i * phi_l) * Y_reg_inv;
            W_MLS_r(k,:) = abs(HR(k,:)) .* exp(1i * phi_r) * Y_reg_inv;
        end
    end
end

% mamnually set the DC bin (use `real()` instead of `abs()`, which causes 
% strong a magnitude errors in the rendering results at low frequencies)
W_MLS_l(1, :) = real(W_MLS_l(2, :));
W_MLS_r(1, :) = real(W_MLS_r(2, :));

% to time domain
if isreal(Y_hor_conj)
    W_MLS_l = [W_MLS_l(1:numPosFreqs, :); flipud(conj(W_MLS_l(2:numPosFreqs-1, :)))];
    W_MLS_r = [W_MLS_r(1:numPosFreqs, :); flipud(conj(W_MLS_r(2:numPosFreqs-1, :)))];
else
    W_MLS_l = getShFreqDomainConjugate(W_MLS_l(1:numPosFreqs, :));
    W_MLS_r = getShFreqDomainConjugate(W_MLS_r(1:numPosFreqs, :));
end
wMlsL = ifft(W_MLS_l);
wMlsR = ifft(W_MLS_r);
if isreal(Y_hor_conj)
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
