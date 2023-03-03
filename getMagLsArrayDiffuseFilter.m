function wAdf = getMagLsArrayDiffuseFilter( ...
    micRadius, micGridAziRad, micGridZenRad, order, fs, len, ...
    shDefinition, shFunction)

if nargin < 8; shFunction = @getSH; end
if nargin < 7 || isempty(shDefinition); shDefinition = 'real'; end

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

% expand to SHs
bn_Hi = sh_repToOrder(bn_Hi.').';

% sample from high order to target grid
Y_Hi_conj = shFunction(simulationOrder, [micGridAziRad, micGridZenRad], shDefinition)';
bn_Lo_Dir = bn_Hi * Y_Hi_conj;

% decompose at low order
Y_Lo = shFunction(order, [micGridAziRad, micGridZenRad], shDefinition);
bn_Lo = bn_Lo_Dir * Y_Lo;

% calculate diffuse-field response
bn_Hi_df = rms(abs(bn_Hi), 2) * sqrt(size(bn_Hi, 2)) / (4*pi);
bn_Lo_df = rms(abs(bn_Lo), 2) * sqrt(size(bn_Lo, 2)) / (4*pi);

% calculate array diffuse-field difference (equivalent to Spherical Head Filter)
W_Alias = bn_Hi_df ./ bn_Lo_df;

% figure('Name', 'ADF RMS abs');
% AKp([ifft(AKsingle2bothSidedSpectrum(bn_Hi_df)), ifft(AKsingle2bothSidedSpectrum(bn_Lo_df))], 'm2d', 'fs', fs);
% AKp(ifft(AKsingle2bothSidedSpectrum(W_Alias)), 'm2d', 'fs', fs, 'c', 'k');
% legend({'DF High', 'DF Low', 'ADF'}, 'Location', 'SouthEast'); drawnow;
% 
% return

% combine with Spherical Head Filter
[~, W_Shf] = getMagLsSphericalHeadFilter(micRadius, order, fs, len);
W_Adf = W_Shf(1:size(W_Alias, 1), :) .* W_Alias;

% transform into time domain
numPosFreqs = nfft/2+1;
W_Adf = [W_Adf(1:numPosFreqs, :); flipud(conj(W_Adf(2:numPosFreqs-1, :)))];
wAdf = ifft(W_Adf);

% shift from zero-phase-like to linear-phase-like
n_shift = nfft/2;
wAdf = applySubsampleDelay(wAdf, n_shift);

% shorten to target length
wAdf = wAdf(n_shift-len/2+1:n_shift+len/2, :);

% fade
fade_win = getFadeWindow(len);
wAdf = wAdf .* fade_win;

% normalize amplitude
max_amp = max(abs(wAdf));
wAdf = wAdf / max_amp;

end
