function sigFiltered = applyRadialFilter(inSig, params)
% apply a radial filter to a SH-domain SMA signal
% Thomas Deppisch, 2021

radFilts = getRadialFilter(params);
% to time domain
radFilts(isnan(radFilts)) = 0;
irRadFilter = ifft([radFilts; conj(flipud(radFilts(2:end-1,:)))], params.nfft, 1, 'symmetric');
irRadFilter = circshift(irRadFilter, params.nfft/2);

%% apply filter
if size(inSig,1) < params.nfft
    disp('applyRadialFilter: short signal, applying zero padding!')
    inSig = [inSig; zeros(params.nfft-size(inSig,1), size(inSig,2))];
end

sigFiltered = fftfilt(sh_repToOrder(irRadFilter.').', inSig);

% remove filter delay
sigFiltered = sigFiltered(params.nfft/2+1:end,:);

    