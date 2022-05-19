function [wMlsL, wMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, micGridZenRad, order, fs, len, applyDiffusenessConst, ...
    shDefinition, shFunction)
% [wMlsL, wMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
%     micRadius, micGridAziRad, micGridZenRad, order, fs, len, applyDiffusenessConst, ...
%     shDefinition, shFunction)
%
% This function returns eMagLS binaural decoding filters.
% For more information about the renderer, please refer to 
% T. Deppisch, H. Helmholz, J. Ahrens, "End-to-End Magnitude Least Squares Binaural Rendering 
% of Spherical Microphone Array Signals," International 3D Audio Conference (I3DA), 2021.
%
% wMlsL                  .. time-domain decoding filter for left ear
% wMlsR                  .. time-domain decoding filter for right ear
% hL                     .. HRIR set for left e ar (numSamples x numDirections)
% hR                     .. HRIR set for right ear (numSamples x numDirections)
% hrirGridAziRad         .. grid azimuth angles in radians of HRIR set (numDirections x 1)
% hrirGridZenRad         .. grid zenith angles in radians of HRIR set (numDirections x 1)
% micRadius              .. radius of SMA
% micGridAziRad          .. SMA grid azimuth angles in radians
% micGridZenRad          .. SMA grid zenith angles in radians
% order                  .. SH output order
% fs                     .. sampling frequency in Hz
% len                    .. desired length of magLS filters
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
% Thomas Deppisch, 2021

if nargin < 13; shFunction = @getSH; end
if nargin < 12 || isempty(shDefinition); shDefinition = 'real'; end

NFFT_MAX_LEN            = 2048; % maxium length of result in samples
SIMULATION_ORDER        = 32; % see `getSMAIRMatrix()`
SIMULATION_WAVE_MODEL   = 'planeWave'; % see `getSMAIRMatrix()`
SIMULATION_ARRAY_TYPE   = 'rigid'; % see `getSMAIRMatrix()`
F_CUT                   = 2000; % transition frequency in Hz
SVD_REGUL_CONST         = 0.01;
REL_FADE_LEN            = 0.15; % relative length of result fading window

assert(len >= size(hL, 1), 'len too short');

numHarmonics = (order+1)^2;
numDirections = size(hL, 2);
fprintf('with @%s("%s") ... ', func2str(shFunction), shDefinition);
Y_Hi_conj = shFunction(SIMULATION_ORDER, [hrirGridAziRad, hrirGridZenRad], shDefinition)';
Y_Lo_pinv = pinv(Y_Hi_conj(1:numHarmonics, :));

nfft = max(2*len, NFFT_MAX_LEN);
f = linspace(0, fs/2, nfft/2+1).';
k_cut = ceil(F_CUT/f(2));
numPosFreqs = length(f);

% zero pad and remove group delay (alternative to applying global phase delay later)
hL(end+1:nfft, :) = 0;
hR(end+1:nfft, :) = 0;
grpDL = grpdelay(hL * real(Y_Lo_pinv(:, 1)), 1, f, fs);
grpDR = grpdelay(hR * real(Y_Lo_pinv(:, 1)), 1, f, fs);
hL = circshift(hL, -round(median(grpDL)));
hR = circshift(hR, -round(median(grpDR)));

% transform into frequency domain
HL = fft(hL);
HR = fft(hR);

% simulate plane wave impinging on SMA 
params.order = order;
params.fs = fs;
params.irLen = nfft;
params.oversamplingFactor = 1;
params.simulateAliasing = true;
params.simulationOrder = SIMULATION_ORDER;
params.radialFilter = 'none';
params.smaRadius = micRadius;
params.smaDesignAziZenRad = [micGridAziRad, micGridZenRad];
params.waveModel = SIMULATION_WAVE_MODEL;
params.arrayType = SIMULATION_ARRAY_TYPE;
params.shDefinition = shDefinition;
params.shFunction = shFunction;
smairMat = getSMAIRMatrix(params);

W_MLS_l = zeros(numPosFreqs, numHarmonics);
W_MLS_r = zeros(numPosFreqs, numHarmonics);
for k = 2:numPosFreqs % leave out first bin, here bn = 0
    pwGrid = smairMat(:,:,k) * Y_Hi_conj;
    [U,S,V] = svd(pwGrid.', 'econ');
    s = diag(S);
    s = 1 ./ max(s, SVD_REGUL_CONST * max(s)); % regularize
    Y_reg_inv = conj(U) * (s .* V.');

    if k >= k_cut % magnitude least-squares
        phiMagLsSmaL = angle(W_MLS_l(k-1,:) * pwGrid);
        phiMagLsSmaR = angle(W_MLS_r(k-1,:) * pwGrid);

        W_MLS_l(k,:) = abs(HL(k,:)) .* exp(1i * phiMagLsSmaL) * Y_reg_inv;
        W_MLS_r(k,:) = abs(HR(k,:)) .* exp(1i * phiMagLsSmaR) * Y_reg_inv;
    else % least-squares
        W_MLS_l(k,:) = HL(k,:) * Y_reg_inv;
        W_MLS_r(k,:) = HR(k,:) * Y_reg_inv;
    end
end

if applyDiffusenessConst 
    % diffuseness constraint after Zaunschirm, Schoerkhuber, Hoeldrich,
    % "Binaural rendering of Ambisonic signals by head-related impulse
    % response time alignment and a diffuseness constraint"
    
    M = zeros(numPosFreqs, 2, 2);
    HCorr = zeros(numPosFreqs, numHarmonics, 2);
    R = zeros(numPosFreqs, 2, 2);
    RHat = zeros(numPosFreqs, 2, 2);
    RCorr = zeros(numPosFreqs, 2, 2);

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

% DC extension
W_MLS_l(1, :) = W_MLS_l(2, :);
W_MLS_r(1, :) = W_MLS_r(2, :);

% fix spectrum (force real against rounding errors)
W_MLS_l(1, :) = real(W_MLS_l(1, :)); % DC bin
W_MLS_r(1, :) = real(W_MLS_r(1, :));
if ~mod(nfft, 2) % is even
    W_MLS_l(end, :) = real(W_MLS_l(end, :)); % Nyquist bin
    W_MLS_r(end, :) = real(W_MLS_r(end, :));
end

% transform into time domain
W_MLS_l = [W_MLS_l; flipud(conj(W_MLS_l(2:end-1, :)))];
W_MLS_r = [W_MLS_r; flipud(conj(W_MLS_r(2:end-1, :)))];
wMlsL = ifft(W_MLS_l);
wMlsR = ifft(W_MLS_r);
if isreal(Y_Hi_conj)
    assert(isreal(wMlsL), 'Resulting decoding filters are not real valued.');
    assert(isreal(wMlsR), 'Resulting decoding filters are not real valued.');
end

% shift from zero-phase-like to linear-phase-like
% and restore initial group-delay difference between ears
n_shift = nfft/2;
wMlsL = circshift(wMlsL, n_shift);
wMlsR = circshift(wMlsR, n_shift);

% shorten to target length
wMlsL = wMlsL(n_shift-len/2+1:n_shift+len/2, :);
wMlsR = wMlsR(n_shift-len/2+1:n_shift+len/2, :);

% fade
n_fadein = round(REL_FADE_LEN * len);
n_fadeout = round(REL_FADE_LEN * len);
hannin = hann(2*n_fadein);
hannout = hann(2*n_fadeout);
fade_win = [hannin(1:end/2); ones(len-(n_fadein+n_fadeout),1); hannout(end/2+1:end)];

wMlsL = wMlsL .* fade_win;
wMlsR = wMlsR .* fade_win;

end
