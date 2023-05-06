function sig = applySubsampleDelay(sig, delay_samples)
    % apply a time delay with sub-sample precision to an input signal
    %
    % This software is licensed under a Non-Commercial Software License
    % (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
    %
    % Hannes Helmholz, 2023

    % generate double-sided delay spectrum
    omega = linspace(0, 0.5, size(sig, 1)/2+1).';
    exp_omega = exp(-1j * 2 * pi * omega .* delay_samples);
    exp_omega(end, :, :) = real(exp_omega(end, :, :)); % fix Nyquist bin
    exp_omega = [exp_omega; flipud(conj(exp_omega(2:end-1, :, :)))];

    % apply complex delay in frequency domain
    Sig = fft(sig);
    Sig = Sig .* exp_omega;
    sig = ifft(Sig);
end
