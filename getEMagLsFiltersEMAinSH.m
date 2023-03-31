function [wMlsL, wMlsR] = getEMagLsFiltersEMAinSH(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, order, fs, len, applyDiffusenessConst, ...
    shDefinition, shFunction)
% [wMlsL, wMlsR] = getEMagLsFiltersEMAinSH(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
%     micRadius, micGridAziRad, order, fs, len, applyDiffusenessConst, ...
%     shDefinition, shFunction)
%
% This function returns eMagLS binaural decoding filters for equatorial microphone arrays.
%
% wMlsL                  .. time-domain decoding filter for left ear
% wMlsR                  .. time-domain decoding filter for right ear
% hL                     .. HRIR set for left e ar (numSamples x numDirections)
% hR                     .. HRIR set for right ear (numSamples x numDirections)
% hrirGridAziRad         .. grid azimuth angles in radians of HRIR set (numDirections x 1)
% hrirGridZenRad         .. grid zenith angles in radians of HRIR set (numDirections x 1)
% micRadius              .. radius of SMA
% micGridAziRad          .. SMA grid azimuth angles in radians
% order                  .. SH output order
% fs                     .. sampling frequency in Hz
% len                    .. desired length of eMagLS filters
% applyDiffusenessConst  .. {true, false}, apply diffuseness constraint,
%                           see Zaunschirm, Schoerkhuber, Hoeldrich,
%                           "Binaural rendering of Ambisonic signals by head-related impulse
%                           response time alignment and a diffuseness constraint"
% shDefinition           .. SH basis type according to utilized shFunction, default: 'real'
% shFunction             .. SH basis function (see testEMagLs.m for example), default: @getSH
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/master/LICENSE for full details).
%
% Hannes Helmholz, 2023

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
f_cut = max(F_CUT_MIN_FREQ, 500 * order); % from N > k
k_cut = ceil(f_cut / f(2));
fprintf('with transition at %d Hz ... ', ceil(f_cut));

fprintf('with @%s("%s") ... ', func2str(shFunction), shDefinition);
% simulate plane wave impinging on EMA
params.returnRawMicSigs = true; % raw mic signals, no SHs!
params.order = order;
params.fs = fs;
params.irLen = nfft;
params.oversamplingFactor = 1;
params.simulateAliasing = true;
params.radialFilter = 'none';
params.smaRadius = micRadius;
params.smaDesignAziZenRad = [micGridAziRad, ones(size(micGridAziRad)) * pi/2];
params.waveModel = SIMULATION_WAVE_MODEL;
params.arrayType = SIMULATION_ARRAY_TYPE;
params.shDefinition = shDefinition;
params.shFunction = shFunction;
emairMat = getSMAIRMatrix(params);
simulationOrder = sqrt(size(emairMat, 2)) - 1;

% sample the microphone responses to PW directions on the EMA from all HRIR
% directions, but mapped to elevation = 0
Y_hor_conj = shFunction(simulationOrder, [hrirGridAziRad, ones(size(hrirGridAziRad)) * pi/2], shDefinition)';
emairDir = pagemtimes(emairMat, Y_hor_conj);
clear emairMat;

% transform into time domain
numPosFreqs = length(f);
emairDir = permute(emairDir, [3, 1, 2]);
emairDir = [emairDir(1:numPosFreqs, :, :); flipud(conj(emairDir(2:numPosFreqs-1, :, :)))];
emairDir_t = ifft(emairDir);
clear emairDir;

% decompose with EMA functions without radial filters
numHarmonics = (order+1)^2;
numDirections = size(hL, 2);
emairDir_sh = zeros(nfft, numHarmonics, numDirections);
for d = 1 : numDirections
    % TODO: Decide how this function should be made available to this repsitory
    % skip radial filtering here (providing 1 as the argument)!!
    emairDir_sh(:, :, d) = get_sound_field_sh_coeffs_from_ema_t( ...
        emairDir_t(:, :, d), 1, order, micGridAziRad.');
end; clear d;
clear emairDir_t;

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
        emairDir_sh(:, :, d) = emairDir_sh(:, :, d) * sh_rot_matrix;
    end
end; clear d euler_matrix sh_rot_matrix;

% transform into freqeuncy domain
pwGridAll = fft(emairDir_sh);
pwGridAll = permute(pwGridAll, [2, 3, 1]);
clear emairDir_sh;

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
for k = 1:numPosFreqs
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
        if k == numPosFreqs && ~mod(nfft, 2) % Nyquist bin, is even
            W_MLS_l(k,:) = real(abs(HL(k,:)) .* exp(1i * phi_l)) * Y_reg_inv;
            W_MLS_r(k,:) = real(abs(HR(k,:)) .* exp(1i * phi_r)) * Y_reg_inv;
        else
            W_MLS_l(k,:) = abs(HL(k,:)) .* exp(1i * phi_l) * Y_reg_inv;
            W_MLS_r(k,:) = abs(HR(k,:)) .* exp(1i * phi_r) * Y_reg_inv;
        end
    end

    if ~isreal(Y_hor_conj) && k > 1 && (k < numPosFreqs || mod(nfft, 2)) % is odd
        warning(strcat('The rendering filters for "complex" SH basis types has not been validated yet.'));

        % negative frequencies below cut in case of complex-valued SHs
        k_neg = nfft-k+2;
        pwGrid = pwGridAll(:,:,k_neg);
        [U, s, V] = svd(pwGrid.', 'econ', 'vector');
        s = 1 ./ max(s, SVD_REGUL_CONST * max(s)); % regularize
        Y_reg_inv = conj(U) * (s .* V.');

        if k < k_cut % least-squares below cut
            W_MLS_l(k_neg,:) = HL(k_neg,:) * Y_reg_inv;
            W_MLS_r(k_neg,:) = HR(k_neg,:) * Y_reg_inv;
        else % magnitude least-squares above cut
            W_MLS_l(k_neg,:) = abs(HL(k_neg,:)) .* exp(1i * -phi_l) * Y_reg_inv;
            W_MLS_r(k_neg,:) = abs(HR(k_neg,:)) .* exp(1i * -phi_r) * Y_reg_inv;
        end
    end
end

if applyDiffusenessConst
    warning('Functionality of the diffuseness constraint has not been verified yet.');

    % diffuseness constraint after Zaunschirm, Schoerkhuber, Hoeldrich,
    % "Binaural rendering of Ambisonic signals by head-related impulse
    % response time alignment and a diffuseness constraint"

    M = zeros(numPosFreqs, 2, 2, 'like', HL);
    HCorr = zeros(numPosFreqs, numHarmonics, 2, 'like', HL);
    R = zeros(numPosFreqs, 2, 2, 'like', HL);
    RHat = zeros(numPosFreqs, 2, 2, 'like', HL);
    RCorr = zeros(numPosFreqs, 2, 2, 'like', HL);

    for ff = 2:numPosFreqs
        % target covariance via original HRTF set
        H = [HL(ff,:); HR(ff,:)];
        R(ff,:,:) = 1/numDirections * (H * H');
        %R(ff,:,:) = H' * W * H;
        R(abs(imag(R)) < 10e-10) = real(R(abs(imag(R)) < 10e-10)); % neglect small imaginary parts
        X = chol(squeeze(R(ff,:,:))); % chol factor of covariance of HRTF set

        % covariance of magLS HRTF set after rendering
        HHat = [W_MLS_l(ff,:); W_MLS_r(ff,:)];
        RHat(ff,:,:) = 1/(4*pi) * (HHat * smairMat(:,:,ff) * smairMat(:,:,ff)' * HHat');
        RHat(abs(imag(RHat)) < 10e-10) = real(RHat(abs(imag(RHat)) < 10e-10));
        XHat = chol(squeeze(RHat(ff,:,:))); % chol factor of magLS HRTF set in SHD

        [U,S,V] = svd(XHat' * X);

        if any(imag(diag(S)) ~= 0) || any(diag(S) < 0)
            warning('negative or complex singular values, pull out negative/complex and factor into left or right singular vector!')
        end

        M(ff,:,:) = V * U' * X / XHat;
        HCorr(ff,:,:) = HHat' * squeeze(M(ff,:,:));

        RCorr(ff,:,:) = 1/(4*pi) * squeeze(HCorr(ff,:,:))' * smairMat(:,:,ff) * smairMat(:,:,ff)' * squeeze(HCorr(ff,:,:));
    end

    W_MLS_l = conj(HCorr(:,:,1));
    W_MLS_r = conj(HCorr(:,:,2));
end

% transform into time domain
if isreal(Y_hor_conj)
    W_MLS_l = [W_MLS_l(1:numPosFreqs, :); flipud(conj(W_MLS_l(2:numPosFreqs-1, :)))];
    W_MLS_r = [W_MLS_r(1:numPosFreqs, :); flipud(conj(W_MLS_r(2:numPosFreqs-1, :)))];
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
