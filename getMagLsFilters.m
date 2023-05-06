function [wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    order, fs, len, applyDiffusenessConst, shDefinition, shFunction)
% [wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
%     order, fs, len, applyDiffusenessConst, shDefinition, shFunction)
%
% This function calculates MagLS binaural decoding filters for head related 
% impulse response data sets.
% For more information, please refer to
%   Schörkhuber, Zaunschirm, and Hoeldrich,
%   “Binaural Rendering of Ambisonic Signals via Magnitude Least Squares,”
%   in Fortschritte der Akustik -- DAGA 2018, 2018, pp. 339–342.
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
% applyDiffusenessConst  .. {true, false}, apply diffuseness constraint, default: false
% shDefinition           .. SH basis type according to utilized shFunction, default: 'real'
% shFunction             .. SH basis function (see testEMagLs.m for example), default: @getSH
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Thomas Deppisch, 2023

if nargin < 10; shFunction = @getSH; end
if nargin < 9 || isempty(shDefinition); shDefinition = 'real'; end
if nargin < 8 || isempty(applyDiffusenessConst); applyDiffusenessConst = false; end

NFFT_MAX_LEN            = 2048; % maxium oversamping length in samples
F_CUT_MIN_FREQ          = 1e3; % minimum transition freqeuncy in Hz
DIFF_CONST_IMAG_THLD    = 1e-9;

% TODO: Implement dealing with HRIRs that are longer than the requested filter
assert(len >= size(hL, 1), 'len too short');

nfft = min(NFFT_MAX_LEN, 2 * len); % apply frequency-domain oversampling
f = linspace(0, fs/2, nfft/2+1).';
numPosFreqs = length(f);
f_cut = max(F_CUT_MIN_FREQ, 500 * order); % from N > k
k_cut = ceil(f_cut / f(2));
fprintf('with transition at %d Hz ... ', ceil(f_cut));

numHarmonics = (order+1)^2;
numDirections = size(hL, 2);
fprintf('with @%s("%s") ... ', func2str(shFunction), shDefinition);
Y_conj = shFunction(order, [hrirGridAziRad, hrirGridZenRad], shDefinition)';
Y_pinv = pinv(Y_conj);

% estimate group delay, zero pad and remove group delay with subsample
% precision (this is an alternative to applying global phase delay later)
grpD = median(cat(3, grpdelay(sum(hL, 2), 1, f, fs), grpdelay(sum(hR, 2), 1, f, fs)));
h = cat(3, hL, hR);
clear hL hR;
h(end+1:nfft, :, :) = 0;
h = applySubsampleDelay(h, -grpD);

w_LS = pagemtimes(h, Y_pinv);

% transform into frequency domain
H = fft(h);
W_MLS = fft(w_LS);

for k = k_cut:numPosFreqs
    phi = angle(pagemtimes(W_MLS(k-1, :, :), Y_conj));

    if k == numPosFreqs && ~mod(nfft, 2) % Nyquist bin, is even
        W_MLS(k, :, :) = pagemtimes(real(abs(H(k, :, :)) .* exp(1i * phi)), Y_pinv);
    else
        % positive frequencies
        W_MLS(k, :, :) = pagemtimes(abs(H(k, :, :)) .* exp(1i * phi), Y_pinv);
        if ~isreal(Y_conj)
            % negative frequencies in case of complex-valued SHs
            k_neg = nfft-k+2;
            W_MLS(k_neg, :, :) = pagemtimes(abs(H(k_neg, :, :)) .* exp(1i * -phi), Y_pinv);
        end
    end
end

if applyDiffusenessConst
    assert(strcmpi(shDefinition, 'real'), ...
        'Diffuseness constraint is not implemented for "%s" SHs yet.', shDefinition);
    % diffuseness constraint after Zaunschirm, Schoerkhuber, Hoeldrich,
    % "Binaural rendering of Ambisonic signals by head-related impulse
    % response time alignment and a diffuseness constraint"

    HCorr = zeros(numPosFreqs, numHarmonics, 2, 'like', H);
    for k = 1:numPosFreqs
        % target covariance via original HRTF set
        HT = squeeze(H(k, :, :));
        RT = 1/numDirections * (HT' * HT);
        RT_small = abs(imag(RT)) < DIFF_CONST_IMAG_THLD;
        RT(RT_small) = real(RT(RT_small)); % neglect small imaginary parts
        XT = chol(RT); % chol factor of covariance of HRTF set

        % covariance of magLS HRTF set
        HHat = squeeze(W_MLS(k, :, :));
        RHat = 1/(4*pi) * (HHat' * HHat);
        RHat_small = abs(imag(RHat)) < DIFF_CONST_IMAG_THLD; % neglect small imaginary parts
        RHat(RHat_small) = real(RHat(RHat_small));
        XHat = chol(RHat); % chol factor of magLS HRTF set in SHD

        [U, s, V] = svd(XHat' * XT, 'econ', 'vector');
        if any(imag(s) ~= 0) || any(s < 0)
            warning('negative or complex singular values, pull out negative/complex and factor into left or right singular vector!')
        end
        M = V * U' * XT / XHat;
        HCorr(k, :, :) = conj(HHat) * M;
    end

    W_MLS = conj(HCorr);
end

% transform into time domain
if isreal(Y_conj)
   W_MLS = [W_MLS(1:numPosFreqs, :, :); flipud(conj(W_MLS(2:numPosFreqs-1, :, :)))];
end
wMls = ifft(W_MLS);
if isreal(Y_conj)
    assert(isreal(wMls), 'Resulting decoding filters are not real valued.');
end

% shift from zero-phase-like to linear-phase-like
% and restore initial group-delay difference between ears
n_shift = nfft/2;
wMls = applySubsampleDelay(wMls, cat(3, n_shift, n_shift + diff(grpD)));

% shorten to target length
wMls = wMls(n_shift-len/2+1:n_shift+len/2, :, :);

% fade
fade_win = getFadeWindow(len);
wMls = wMls .* fade_win;

wMlsL = wMls(:, :, 1);
wMlsR = wMls(:, :, 2);

end
