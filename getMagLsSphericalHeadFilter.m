function [wShf, W_Shf] = getMagLsSphericalHeadFilter(micRadius, order, fs, len)

NFFT_MAX_LEN            = 2048; % maxium oversamping length in samples
SIMULATION_ARRAY_TYPE   = 'rigid'; % see `sphModalCoeffs()`
SIMULATION_MIC_DIR      = 0; % see `sphModalCoeffs()`
SIMULATION_C            = 343; % speed of sound in m/s

nfft = min(2*len, NFFT_MAX_LEN); % apply frequency-domain oversampling
f = linspace(0, fs/2, nfft/2+1).';
kr = 2*pi*f/SIMULATION_C * micRadius;
simulationOrder = ceil(fs * pi * micRadius / SIMULATION_C);

% get modal radial filters
bn_Hi = sphModalCoeffs(simulationOrder, kr, SIMULATION_ARRAY_TYPE, SIMULATION_MIC_DIR);
bn_Lo = bn_Hi(:, 1:order+1);

% expand to SHs
bn_Hi = sh_repToOrder(bn_Hi.').';
bn_Lo = sh_repToOrder(bn_Lo.').';

% calculate diffuse-field response
bn_Hi_df = rms(abs(bn_Hi), 2) * sqrt(size(bn_Hi, 2)) / (4*pi);
bn_Lo_df = rms(abs(bn_Lo), 2) * sqrt(size(bn_Lo, 2)) / (4*pi);

% calculate scattering diffuse-field difference (Spherical Head Filter)
W_Shf = bn_Hi_df ./ bn_Lo_df;

% figure('Name', 'SHF RMS abs');
% AKp([ifft(AKsingle2bothSidedSpectrum(bn_Hi_df)), ifft(AKsingle2bothSidedSpectrum(bn_Lo_df))], 'm2d', 'fs', fs);
% AKp(ifft(AKsingle2bothSidedSpectrum(W_Shf)), 'm2d', 'fs', fs, 'c', 'k');
% legend({'DF High', 'DF Low', 'SHF'}, 'Location', 'SouthEast'); drawnow;
% 
% return

% invert spectrum
W_Shf = 1 ./ W_Shf;

% transform into time domain
numPosFreqs = nfft/2+1;
W_Shf = [W_Shf(1:numPosFreqs, :); flipud(conj(W_Shf(2:numPosFreqs-1, :)))];
wShf = ifft(W_Shf);

% shift from zero-phase-like to linear-phase-like
n_shift = nfft/2;
wShf = applySubsampleDelay(wShf, n_shift);

% shorten to target length
wShf = wShf(n_shift-len/2+1:n_shift+len/2, :);

% fade
fade_win = getFadeWindow(len);
wShf = wShf .* fade_win;

% normalize amplitude
max_amp = max(abs(wShf));
wShf = wShf / max_amp;

end
