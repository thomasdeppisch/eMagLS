function [wMlsL, wMlsR] = getEMagLs2Filters(hL, hR, hrirGridAziRad, hrirGridZenRad, micRadius, micGridAziRad, micGridZenRad, fs, len, applyDiffusenessConst)
% [wMlsL, wMlsR] = getEMagLs2Filters(hL, hR, hrirGridAziRad, hrirGridZenRad, micRadius, micGridAziRad, micGridZenRad, fs, len, applyDiffusenessConst)
%
% This function returns eMagLS2 binaural decoding filters.
% For more information about the renderer, please refer to 
% T. Deppisch, H. Helmholz, J. Ahrens, "End-to-End Magnitude Least Squares Binaural Rendering 
% of Spherical Microphone Array Signals," International 3D Audio Conference (I3DA), 2021.
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
% fs                     .. sampling frequency in Hz
% len                    .. desired length of magLS filters
% applyDiffusenessConst  .. {true, false}, apply diffuseness constraint,
%                           see Zaunschirm, Schoerkhuber, Hoeldrich,
%                           "Binaural rendering of Ambisonic signals by head-related impulse
%                           response time alignment and a diffuseness constraint"
%
% This software is licensed under a Non-Commercial Software License 
% (see https://github.com/thomasdeppisch/eMagLS/blob/master/LICENSE for full details).
%
% Thomas Deppisch, 2021

shDefinition = 'complex'; % real or complex

if (len < size(hL,1))
    error('len too short')
end

numMics = length(micGridAziRad);
numDirections = size(hL,2);
nfft = max(2048,2*len);
simulationOrder = 32;
YHi = getSH(simulationOrder, [hrirGridAziRad,hrirGridZenRad], shDefinition).';

f_cut = 2000; % transition frequency 
f = linspace(0,fs/2,nfft/2+1);
k_cut = round(f_cut/f(2) + 1);

% zero pad and remove group delay (alternative to applying global phase delay later)
grpDL = grpdelay(mean(hL,2), 1, f, fs);
grpDR = grpdelay(mean(hR,2), 1, f, fs);

hL = circshift([hL; zeros(nfft - size(hL, 1), size(hL, 2))], -round(median(grpDL)));
hR = circshift([hR; zeros(nfft - size(hR, 1), size(hR, 2))], -round(median(grpDR)));

HL = fft(hL,nfft);
HR = fft(hR,nfft);

numPosFreqs = nfft/2+1;

% simulate plane wave impinging on SMA 
params.returnRawMicSigs = true; % raw mic signals, no SHs!
params.fs = fs;
params.irLen = nfft;
params.oversamplingFactor = 1;
params.simulateAliasing = true;
params.simulationOrder = simulationOrder;
params.radialFilter = 'none';
params.smaRadius = micRadius;
params.smaDesignAziZenRad = [micGridAziRad, micGridZenRad];
params.waveModel = 'planeWave';
params.arrayType = 'rigid';
params.shDefinition = shDefinition;
smairMat = getSMAIRMatrix(params);

svdRegulConst = 0.01;

W_LS_l = zeros(nfft, numMics);
W_LS_r = zeros(nfft, numMics);
W_MLS_l = zeros(nfft, numMics);
W_MLS_r = zeros(nfft, numMics);
for k = 1:numPosFreqs % leave out first bin, here bn = 0
    pwGrid = smairMat(:,:,k) * conj(YHi);
    [U,S,V] = svd(pwGrid.','econ');
    s = diag(S);
    s = 1 ./ max(s, svdRegulConst * max(s)); % regularize
    regInvY = (V .* s.') * U';
    
    W_LS_l(k,:) = (regInvY * HL(k,:).').';
    W_LS_r(k,:) = (regInvY * HR(k,:).').';

    % calculate negative frequencies (important in case of complex-valued
    % SHs)
    if k > 1 && k < numPosFreqs
        pwGridNeg = smairMat(:,:,end-k+2) * conj(YHi);
        [UNeg,SNeg,VNeg] = svd(pwGridNeg.','econ');
        sNeg = diag(SNeg);
        sNeg = 1 ./ max(sNeg, svdRegulConst * max(sNeg)); % regularize
        regInvYNeg = (VNeg .* sNeg.') * UNeg';

        W_LS_l(end-k+2,:) = (regInvYNeg * conj(HL(k,:)).').';
        W_LS_r(end-k+2,:) = (regInvYNeg * conj(HR(k,:)).').';
    end
    
    if k >= k_cut
        phiMagLsSmaL = angle(pwGrid.' * W_MLS_l(k-1,:).');
        phiMagLsSmaR = angle(pwGrid.' * W_MLS_r(k-1,:).');
        
        W_MLS_l(k,:) = regInvY * (abs(HL(k,:)).' .* exp(1i * phiMagLsSmaL));
        W_MLS_r(k,:) = regInvY * (abs(HR(k,:)).' .* exp(1i * phiMagLsSmaR));

        if k > 1 && k < numPosFreqs
            phiMagLsSmaLNeg = angle(pwGridNeg.' * conj(W_MLS_l(k-1,:)).');
            phiMagLsSmaRNeg = angle(pwGridNeg.' * conj(W_MLS_r(k-1,:)).');
            
            W_MLS_l(end-k+2,:) = regInvYNeg * (abs(HL(k,:)).' .* exp(1i * phiMagLsSmaLNeg));
            W_MLS_r(end-k+2,:) = regInvYNeg * (abs(HR(k,:)).' .* exp(1i * phiMagLsSmaRNeg));
        end

        % Nyquist
        if k == numPosFreqs
            W_MLS_l(k,:) = regInvY * real(abs(HL(k,:)).' .* exp(1i * phiMagLsSmaL));
            W_MLS_r(k,:) = regInvY * real(abs(HR(k,:)).' .* exp(1i * phiMagLsSmaR));
        end
    else
        W_MLS_l(k,:) = W_LS_l(k,:);
        W_MLS_r(k,:) = W_LS_r(k,:);

        if k > 1 && k < numPosFreqs
            W_MLS_l(end-k+2,:) = W_LS_l(end-k+2,:);
            W_MLS_r(end-k+2,:) = W_LS_r(end-k+2,:);
        end
    end
end

if applyDiffusenessConst
    % diffuseness constraint after Zaunschirm, Schoerkhuber, Hoeldrich,
    % "Binaural rendering of Ambisonic signals by head-related impulse
    % response time alignment and a diffuseness constraint"
    
    M = zeros(numPosFreqs, 2, 2);
    HCorr = zeros(numPosFreqs, numMics, 2);
    R = zeros(numPosFreqs, 2, 2);
    RHat = zeros(numPosFreqs, 2, 2);
    RCorr = zeros(numPosFreqs, 2, 2);

    for ff = 2:numPosFreqs
        % target covariance via original HRTF set
        H = [HL(ff,:); HR(ff,:)];
        R(ff,:,:) = 1/numDirections * (H * H');
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

W_MLS_l(1,:) = real(W_MLS_l(2,:)); % DC extension
W_MLS_r(1,:) = real(W_MLS_r(2,:));
W_MLS_l(end,:) = real(W_MLS_l(end,:));
W_MLS_r(end,:) = real(W_MLS_r(end,:));
W_MLS_l = [W_MLS_l; flipud(conj(W_MLS_l(2:end-1,:)))];
wMlsL = ifft(W_MLS_l,nfft);
W_MLS_r = [W_MLS_r; flipud(conj(W_MLS_r(2:end-1,:)))];
wMlsR = ifft(W_MLS_r,nfft);

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

