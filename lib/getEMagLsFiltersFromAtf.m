function [wMlsL, wMlsR] = getEMagLsFiltersFromAtf(hL, hR, hrirGridAziZenRad, atfIrs, atfGridAziZenRad, fs, filterLen, fTrans)
% [wMlsL, wMlsR] = getEMagLsFiltersFromAtf(hL, hR, hrirGridAziZenRad, atfIrs, atfGridAziZenRad, fs, filterLen, fTrans)
%
% This function calculates eMagLS2 binaural decoding filters for arbitrary
% microphone arrays based on array transfer functions (ATFs).
%
% For more information, please refer to
%   Thomas Deppisch, Nils Meyer-Kahlen, Sebastia Amengual Gari, 
%   "Blind Identification of Binaural Room Impulse Responses from Smart Glasses", 
%   IEEE/ACM Transactions on Audio, Speech, and Language Processing, vol. 32, pp. 4052-4065, 2024, 
%   doi: 10.1109/TASLP.2024.3454964.
%
% wMlsL                  .. time-domain decoding filter for left ear
% wMlsR                  .. time-domain decoding filter for right ear
% hL                     .. HRIR set for left ear (numSamples x numDirections)
% hR                     .. HRIR set for right ear (numSamples x numDirections)
% hrirGridAziZenRad      .. grid azimuth and zenith angles in radians of HRIR set (numDirections x 2)
% atfIrs                 .. set of array transfer functions (ATFs)
% atfGridAziZenRad       .. grid azimuth and zenith angles in radians of ATFs
% fs                     .. sampling frequency in Hz
% filterLen              .. desired length of eMagLS2 filters
% fTrans                 .. frequency above which the magnitude will be optimized
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Thomas Deppisch, 2023

NFFT_MAX_LEN            = 2048; % maximum oversampling length in samples
SVD_REGUL_CONST         = 0.01;

% TODO: Implement dealing with HRIRs that are longer than the requested filter
assert(filterLen >= size(hL, 1), 'len too short');

nfft = min(NFFT_MAX_LEN, 2 * filterLen); % apply frequency-domain oversampling
f = linspace(0, fs/2, nfft/2+1).';
numPosFreqs = length(f);
kTrans = ceil(fTrans / f(2));

numMics = size(atfIrs,2);

% zero pad and remove group delay
hL(end+1:nfft, :) = 0;
hR(end+1:nfft, :) = 0;
grpDL = median(grpdelay(sum(hL, 2), 1, f, fs));
grpDR = median(grpdelay(sum(hR, 2), 1, f, fs));

hL = circshift([hL; zeros(nfft - size(hL, 1), size(hL, 2))], -round(grpDL));
hR = circshift([hR; zeros(nfft - size(hR, 1), size(hR, 2))], -round(grpDR));

% transform into frequency domain
HL = fft(hL, nfft);
HR = fft(hR, nfft);
atfs = fft(atfIrs, nfft);

[hrirGridCart(:,1), hrirGridCart(:,2), hrirGridCart(:,3)] = sph2cart(hrirGridAziZenRad(:,1),pi/2-hrirGridAziZenRad(:,2),1);
[atfGridCart(:,1), atfGridCart(:,2), atfGridCart(:,3)] = sph2cart(atfGridAziZenRad(:,1),pi/2-atfGridAziZenRad(:,2),1);

% match the HRTF grid and the ATF grid:
% take all points from the smaller grids and find closest directions from
% the larger grid
[numDirections, idxOfSmallerGrid] = min([size(hL, 2), size(atfIrs,3)]);
if idxOfSmallerGrid == 1
    % HRTF grid is smaller
    dirGridAziZenRad = hrirGridAziZenRad;
    dirGridCart = hrirGridCart;
    gridToMatchCart = atfGridCart;
    HLMatched = HL;
    HRMatched = HR;
    atfsMatched = zeros(numPosFreqs, numMics, numDirections);
else
    % ATF grid is smaller
    dirGridAziZenRad = atfGridAziZenRad;
    dirGridCart = atfGridCart;
    gridToMatchCart = hrirGridCart;
    HLMatched = zeros(numPosFreqs, numDirections);
    HRMatched = zeros(numPosFreqs, numDirections);
    atfsMatched = atfs;
end

gridPointAngDistDeg = zeros(numDirections,1);
for ii = 1:numDirections
    % match
    [~, closestDirIdx] = min(sqrt(sum((gridToMatchCart - dirGridCart(ii,:)).^2, 2))); % find closes direction
    gridPointAngDistDeg(ii) = acos(dirGridCart(ii,:) * gridToMatchCart(closestDirIdx,:)') * 180/pi;

    if idxOfSmallerGrid == 1
        % HRTF grid is smaller
        atfsMatched(:, :, ii) = atfs(1:numPosFreqs, :, closestDirIdx);
    else
        % ATF grid is smaller
        HLMatched(:, ii) = HL(1:numPosFreqs, closestDirIdx);
        HRMatched(:, ii) = HR(1:numPosFreqs, closestDirIdx);
    end
end
disp(['Matching HRTF and ATF grids, average grid deviation: ' num2str(mean(gridPointAngDistDeg)) ' deg'])

W_MLS_l = zeros(numPosFreqs, numMics, 'like', HL);
W_MLS_r = zeros(numPosFreqs, numMics, 'like', HL);
for k = 2:numPosFreqs
    pwGrid = squeeze(atfsMatched(k, :, :)); %  numMics x numDirections
    [U, s, V] = svd(pwGrid.', 'econ', 'vector');
    s = 1 ./ max(s, SVD_REGUL_CONST * max(s)); % regularize
    Y_reg_inv = conj(U) * (s .* V.'); 

    if k < kTrans % least-squares below cut
        W_MLS_l(k,:) = HLMatched(k,:) * Y_reg_inv;
        W_MLS_r(k,:) = HRMatched(k,:) * Y_reg_inv;
    else % magnitude least-squares above cut
        phi_l = angle(W_MLS_l(k-1,:) * pwGrid);
        phi_r = angle(W_MLS_r(k-1,:) * pwGrid);
        if k == numPosFreqs && ~mod(nfft, 2) % Nyquist bin, is even
            W_MLS_l(k,:) = real(abs(HLMatched(k,:)) .* exp(1i * phi_l)) * Y_reg_inv;
            W_MLS_r(k,:) = real(abs(HRMatched(k,:)) .* exp(1i * phi_r)) * Y_reg_inv;
        else
            W_MLS_l(k,:) = abs(HLMatched(k,:)) .* exp(1i * phi_l) * Y_reg_inv;
            W_MLS_r(k,:) = abs(HRMatched(k,:)) .* exp(1i * phi_r) * Y_reg_inv;
        end
    end
end


% mamnually set the DC bin (use `real()` instead of `abs()`, which causes
% strong a magnitude errors in the rendering results at low frequencies)
W_MLS_l(1, :) = real(W_MLS_l(2, :));
W_MLS_r(1, :) = real(W_MLS_r(2, :));

% transform into time domain
W_MLS_l = [W_MLS_l(1:numPosFreqs, :); flipud(conj(W_MLS_l(2:numPosFreqs-1, :)))];
W_MLS_r = [W_MLS_r(1:numPosFreqs, :); flipud(conj(W_MLS_r(2:numPosFreqs-1, :)))];
wMlsL = ifft(W_MLS_l);
wMlsR = ifft(W_MLS_r);

% shift from zero-phase-like to linear-phase-like
% and restore initial group-delay difference between ears
n_shift = round(nfft/2);
wMlsL = circshift(wMlsL, n_shift);
wMlsR = circshift(wMlsR, n_shift);

% shorten to target length
wMlsL = wMlsL(n_shift-filterLen/2+1:n_shift+filterLen/2, :);
wMlsR = wMlsR(n_shift-filterLen/2+1:n_shift+filterLen/2, :);

% fade
n_fadein = round(0.15 * filterLen);
n_fadeout = round(0.15 * filterLen);
hannin = hann(2*n_fadein);
hannout = hann(2*n_fadeout);
fade_win = [hannin(1:end/2); ones(filterLen-(n_fadein+n_fadeout),1); hannout(end/2+1:end)];
wMlsL = wMlsL .* fade_win;
wMlsR = wMlsR .* fade_win;

end