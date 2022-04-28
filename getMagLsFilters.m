function [wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    order, fs, len, applyDiffusenessConst, shDefinition)
% [wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
%     order, fs, len, applyDiffusenessConst, shDefinition)
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
% shDefinition           .. {'real', 'complex'}, SH basis type, default: 'real'
%
% This software is licensed under a Non-Commercial Software License 
% (see https://github.com/thomasdeppisch/eMagLS/blob/master/LICENSE for full details).
%
% Thomas Deppisch, 2021

if nargin < 9; shDefinition = 'real'; end

if (len < size(hL,1))
    error('len too short')
end

nfft = max(2*len,2048);
Y = getSH(order, [hrirGridAziRad,hrirGridZenRad], shDefinition);
pinvY = pinv(Y);

f_cut = 2000; % transition frequency 
f = linspace(0,fs/2,nfft/2+1);
k_cut = round(f_cut/f(2) + 1);

% zero pad and remove delay (alternative to applying global phase delay later)
grpDL = grpdelay(hL * pinvY(1,:)', 1, f, fs);
grpDR = grpdelay(hR * pinvY(1,:)', 1, f, fs);

hL = circshift([hL; zeros(nfft - size(hL, 1), size(hL, 2))], -round(median(grpDL)));
hR = circshift([hR; zeros(nfft - size(hR, 1), size(hR, 2))], -round(median(grpDR)));
HL = fft(hL,nfft);
HR = fft(hR,nfft);

w_LS_l = hL * pinvY';
w_LS_r = hR * pinvY';

W_LS_l = fft(w_LS_l,nfft);
W_LS_r = fft(w_LS_r,nfft);

numPosFreqs = nfft/2+1;

W_MLS_l = W_LS_l(1:numPosFreqs,:);
W_MLS_r = W_LS_r(1:numPosFreqs,:);
for k = k_cut:numPosFreqs
    phi_l = angle(W_MLS_l(k-1,:) * Y');
    W_MLS_l(k,:) = (abs(HL(k,:)) .* exp(1i * phi_l)) * pinvY';
    
    phi_r = angle(W_MLS_r(k-1,:) * Y');
    W_MLS_r(k,:) = (abs(HR(k,:)) .* exp(1i * phi_r)) * pinvY';
end

if applyDiffusenessConst
    % diffuseness constraint after Zaunschirm, Schoerkhuber, Hoeldrich,
    % "Binaural rendering of Ambisonic signals by head-related impulse
    % response time alignment and a diffuseness constraint"
    
    numDirections = size(hL,2);
    numHarmonics = (order+1)^2;
    
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

W_MLS_l = [W_MLS_l; flipud(conj(W_MLS_l(2:end-1,:)))];
wMlsL = ifft(W_MLS_l,nfft,'symmetric');
W_MLS_r = [W_MLS_r; flipud(conj(W_MLS_r(2:end-1,:)))];
wMlsR = ifft(W_MLS_r,nfft,'symmetric');

% shorten, shift
n_shift = nfft/2;
wMlsL = circshift(wMlsL, n_shift);
wMlsR = circshift(wMlsR, n_shift);
wMlsL = wMlsL(nfft/2-len/2+1:nfft/2+len/2,:);
wMlsR = wMlsR(nfft/2-len/2+1:nfft/2+len/2,:);

% fade
n_fadein = round(0.15 * len);
n_fadeout = round(0.15 * len);
hannin = hann(2*n_fadein);
hannout = hann(2*n_fadeout);
fade_win = [hannin(1:end/2); ones(len-(n_fadein+n_fadeout),1); hannout(end/2+1:end)];

wMlsL = wMlsL .* fade_win;
wMlsR = wMlsR .* fade_win;

end
