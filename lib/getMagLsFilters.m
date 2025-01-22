function [wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
                                          order, fs, len, shDefinition, shFunction)
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
% shDefinition           .. SH basis type according to utilized shFunction, default: 'real'
% shFunction             .. SH basis function (see testEMagLs.m for example), default: @getSH
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Thomas Deppisch, 2023

if nargin < 9; shFunction = @getSH; end
if nargin < 8 || isempty(shDefinition); shDefinition = 'real'; end

NFFT_MAX_LEN            = 2048; % maxium oversamping length in samples
F_CUT_MIN_FREQ          = 1e3; % minimum transition freqeuncy in Hz

% TODO: Implement dealing with HRIRs that are longer than the requested filter
assert(len >= size(hL, 1), 'HRIR len too short');

nfft = min(NFFT_MAX_LEN, 2 * len); % apply frequency-domain oversampling
f = linspace(0, fs/2, nfft/2+1).'; % assuming an even fft length
numPosFreqs = length(f);
f_cut = max(F_CUT_MIN_FREQ, 500 * order); % from N > k
k_cut = ceil(f_cut / f(2));

numHarmonics = (order+1)^2;
numDirections = size(hL, 2);
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

% to frequency domain
H = fft(h,nfft);
W_MLS = fft(w_LS,nfft);

for k = k_cut:numPosFreqs
    phi = angle(pagemtimes(W_MLS(k-1, :, :), Y_conj));

    if k == numPosFreqs % Nyquist bin, is even
        W_MLS(k, :, :) = pagemtimes(real(abs(H(k, :, :)) .* exp(1i * phi)), Y_pinv);
    else
        W_MLS(k, :, :) = pagemtimes(abs(H(k, :, :)) .* exp(1i * phi), Y_pinv);
    end
end

% to time domain
if isreal(Y_conj)
   W_MLS = [W_MLS(1:numPosFreqs, :, :); flipud(conj(W_MLS(2:numPosFreqs-1, :, :)))];
else
   W_MLS(:,:,1) = getShFreqDomainConjugate(W_MLS(1:numPosFreqs,:,1));
   W_MLS(:,:,2) = getShFreqDomainConjugate(W_MLS(1:numPosFreqs,:,2));
end
wMls = ifft(W_MLS);

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
