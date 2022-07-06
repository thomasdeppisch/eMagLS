function sigFiltered = applyRadialFilter(inSig, params)
% apply a radial filter to a SH-domain SMA signal
% Thomas Deppisch, 2021

radFilts = getRadialFilter(params);
% to time domain
radFilts(isnan(radFilts)) = 0;
irRadFilter = ifft([radFilts; conj(flipud(radFilts(2:end-1,:)))], 'symmetric');

radFiltLen = size(irRadFilter,1);
irRadFilter = circshift(irRadFilter, radFiltLen/2);

% apply fade
fadeLenSmp = 20;
fadeWin = hann(2*fadeLenSmp);
irRadFilter(1:fadeLenSmp,:) = irRadFilter(1:fadeLenSmp,:) .* fadeWin(1:fadeLenSmp);
irRadFilter(end-fadeLenSmp+1:end,:) = irRadFilter(end-fadeLenSmp+1:end,:) .* fadeWin(end-fadeLenSmp+1:end,:);

%% apply filter
if size(inSig,1) < radFiltLen
    disp('applyRadialFilter: short signal, applying zero padding!')
    inSig = [inSig; zeros(params.nfft-size(inSig,1), size(inSig,2))];
end

sigFiltered = fftfilt(sh_repToOrder(irRadFilter.').', inSig);

% remove filter delay
sigFiltered = sigFiltered(radFiltLen/2+1:end,:);

    