function [wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    order, fs, len, applyDiffusenessConst, shDefinition, shFunction)
% [wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
%     order, fs, len, applyDiffusenessConst, shDefinition, shFunction)
%
% calculates magLS binaural decoding filters
% see Schoerkhuber, Zaunschirm, Hoeldrich,
% "Binaural Rendering of Ambisonic Signals via Magnitude Least Squares"
%
% wMlsL                  .. time-domain decoding filter for left ear
% wMlsR                  .. time-domain decoding filter for right ear
% hL                     .. HRIR set for left ear (numSamples x numDirections)
% hR                     .. HRIR set for right ear (numSamples x numDirections)
% hrirGridAziRad         .. grid azimuth angles in radians of HRIR set (numDirections x 1)
% hrirGridZenRad         .. grid zenith angles in radians of HRIR set (numDirections x 1)
% order                  .. SH output order
% fs                     .. sampling frequency in Hz
% len                    .. desired length of magLS filters
% applyDiffuseFieldConst .. {true, false}, apply diffuse-field constraint,
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

if nargin < 10; shFunction = @getSH; end
if nargin < 9 || isempty(shDefinition); shDefinition = 'real'; end

NFFT_MAX_LEN    = 2048; % maxium length of result in samples
REL_FADE_LEN    = 0.15; % relative length of result fading window

assert(len >= size(hL, 1), 'len too short');

nfft = max(2*len, NFFT_MAX_LEN);
f = linspace(0, fs/2, nfft/2+1).';
numPosFreqs = length(f);
f_cut = 500 * order; % from N > k
k_cut = ceil(f_cut / f(2));
fprintf('with transition at %d Hz ... ', ceil(f_cut));

numHarmonics = (order+1)^2;
numDirections = size(hL, 2);
fprintf('with @%s("%s") ... ', func2str(shFunction), shDefinition);
Y_conj = shFunction(order, [hrirGridAziRad, hrirGridZenRad], shDefinition)';
Y_pinv = pinv(Y_conj);

% zero pad and remove group delay with subsample precision
% (alternative to applying global phase delay later)
hL(end+1:nfft, :) = 0;
hR(end+1:nfft, :) = 0;
grpDL = median(grpdelay(hL * Y_pinv(:, 1), 1, f, fs));
grpDR = median(grpdelay(hR * Y_pinv(:, 1), 1, f, fs));
hL = applySubsampleDelay(hL, -grpDL);
hR = applySubsampleDelay(hR, -grpDR);

w_LS_l = hL * Y_pinv;
w_LS_r = hR * Y_pinv;

% transform into frequency domain
HL = fft(hL);
HR = fft(hR);
W_LS_l = fft(w_LS_l);
W_LS_r = fft(w_LS_r);

W_MLS_l = W_LS_l;
W_MLS_r = W_LS_r;
for k = k_cut:numPosFreqs
    phi_l = angle(W_MLS_l(k-1,:) * Y_conj);
    phi_r = angle(W_MLS_r(k-1,:) * Y_conj);

    if k == numPosFreqs && ~mod(nfft, 2) % Nyquist bin, is even
        W_MLS_l(k,:) = real(abs(HL(k,:)) .* exp(1i * phi_l)) * Y_pinv;
        W_MLS_r(k,:) = real(abs(HR(k,:)) .* exp(1i * phi_r)) * Y_pinv;
    else
        % positive frequencies
        W_MLS_l(k,:) = abs(HL(k,:)) .* exp(1i * phi_l) * Y_pinv;
        W_MLS_r(k,:) = abs(HR(k,:)) .* exp(1i * phi_r) * Y_pinv;
        if ~isreal(Y_conj)
            % negative frequencies in case of complex-valued SHs
            k_neg = nfft-k+2;
            W_MLS_l(k_neg,:) = abs(HL(k_neg,:)) .* exp(1i * -phi_l) * Y_pinv;
            W_MLS_r(k_neg,:) = abs(HR(k_neg,:)) .* exp(1i * -phi_r) * Y_pinv;
        end
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

    for ff = 1:numPosFreqs
        % target covariance via original HRTF set
        H = [HL(ff,:); HR(ff,:)];
        R(ff,:,:) = 1/numDirections * (H * H');
        R(abs(imag(R)) < 10e-10) = real(R(abs(imag(R)) < 10e-10)); % neglect small imaginary parts
        X = chol(squeeze(R(ff,:,:))); % chol factor of covariance of HRTF set

        % covariance of magLS HRTF set
        HHat = [W_MLS_l(ff,:); W_MLS_r(ff,:)];
        RHat(ff,:,:) = 1/(4*pi) * (HHat * HHat');
        RHat(abs(imag(RHat)) < 10e-10) = real(RHat(abs(imag(RHat)) < 10e-10));
        XHat = chol(squeeze(RHat(ff,:,:))); % chol factor of magLS HRTF set in SHD

        [U,S,V] = svd(XHat' * X);

        if any(imag(diag(S)) ~= 0) || any(diag(S) < 0)
            warning('negative or complex singular values, pull out negative/complex and factor into left or right singular vector!')
        end

        M(ff,:,:) = V * U' * X / XHat;
        HCorr(ff,:,:) = HHat' * squeeze(M(ff,:,:));

        RCorr(ff,:,:) = 1/(4*pi) * squeeze(HCorr(ff,:,:))' * squeeze(HCorr(ff,:,:));
    end
    
    W_MLS_l = conj(HCorr(:,:,1));
    W_MLS_r = conj(HCorr(:,:,2));
end

% transform into time domain
if isreal(Y_conj)
   W_MLS_l = [W_MLS_l(1:numPosFreqs, :); flipud(conj(W_MLS_l(2:numPosFreqs-1, :)))];
   W_MLS_r = [W_MLS_r(1:numPosFreqs, :); flipud(conj(W_MLS_r(2:numPosFreqs-1, :)))];
end
wMlsL = ifft(W_MLS_l);
wMlsR = ifft(W_MLS_r);
if isreal(Y_conj)
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
n_fadein = round(REL_FADE_LEN * len);
n_fadeout = round(REL_FADE_LEN * len);
hannin = hann(2*n_fadein);
hannout = hann(2*n_fadeout);
fade_win = [hannin(1:end/2); ones(len-(n_fadein+n_fadeout),1); hannout(end/2+1:end)];

wMlsL = wMlsL .* fade_win;
wMlsR = wMlsR .* fade_win;

end
