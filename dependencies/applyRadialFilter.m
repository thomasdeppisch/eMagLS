function sigFiltered = applyRadialFilter(inSig, params)
% apply a radial filter to a SH-domain SMA signal
% Thomas Deppisch, 2021

radFilts = getRadialFilter(params);
radFilts(isnan(radFilts)) = 0;

% transform into time domain
irRadFilter = ifft([radFilts; flipud(conj(radFilts(2:end-1,:)))]);

% make causal
irRadFilter = applySubsampleDelay(irRadFilter, params.nfft/2);

% fade
fade_win = getFadeWindow(params.nfft, 0.05);
irRadFilter = irRadFilter .* fade_win;

%% apply filter
if size(inSig,1) < params.nfft
    disp('applyRadialFilter: short signal, applying zero padding!');
    inSig(end+1:params.nfft, :) = 0;
end

sigFiltered = fftfilt(sh_repToOrder(irRadFilter.').', inSig);

% remove filter delay
sigFiltered = sigFiltered(params.nfft/2+1:end,:);

end
