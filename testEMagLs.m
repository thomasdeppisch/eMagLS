% This script demonstrates how the different eMagLS renderers can be used
% to render recordings from spherical microphone arrays (SMAs) or
% equatorial microphone arrays (EMAs).
%
% For more information, please refer to
%   T. Deppisch, H. Helmholz, and J. Ahrens,
%   “End-to-End Magnitude Least Squares Binaural Rendering of Spherical Microphone Array Signals,”
%   in 2021 Immersive and 3D Audio: from Architecture to Automotive (I3DA), 2021, pp. 1–7. doi: 10.1109/I3DA48870.2021.9610864.
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Thomas Deppisch, 2025

clear; clc; close all;

addpath(genpath('dependencies/'));
addpath(genpath('lib/'));

%% configuration
filterLen = 512; % in samples
shDefinition = 'real';
order = 4;

[hrirFile, hrirUrl] = deal('./resources/HRIR_L2702.mat', ...
    'https://zenodo.org/record/3928297/files/HRIR_L2702.mat');

% load simulated RIRs, captured with an SMA and an EMA
srirSmaStruct = load('./resources/rirSimSma_8cm_32ch_rigid_8x6x4m_278ms.mat');
srirEmaStruct = load('./resources/rirSimEma_8cm_60ch_rigid_8x6x4m_278ms.mat');
sigPath = './resources/decemberTour.wav';

smaGridAziRad = srirSmaStruct.micsAziZenRad(:,1);
smaGridZenRad = srirSmaStruct.micsAziZenRad(:,2);
smaRadius = srirSmaStruct.smaRadius;

emaGridAziRad = srirEmaStruct.micsAziZenRad(:,1);
emaRadius = srirEmaStruct.smaRadius;

%% load data
fprintf('Downloading HRIR dataset ... \n');
if isfile(hrirFile)
    fprintf('already exists ... skipped.\n');
else
    downloadAndExtractFile(hrirFile, hrirUrl);
end

fprintf('Loading file "%s" ... \n', hrirFile);
% the data fields will be incorrect if a different HRIR data set is used
load(hrirFile);
if isa(HRIR_L2702, 'uint32') || isempty(HRIR_L2702)
    % MIRO class does not exist on the system and will be downloaded
    fprintf('Downloading MIRO class ... \n');
    downloadAndExtractFile('dependencies/miro.m', 'https://zenodo.org/record/3928297/files/miro.m');
    % Retry loading file
    fprintf('Loading file "%s" ... \n', hrirFile);
    load(hrirFile);
end
hL = double(HRIR_L2702.irChOne);
hR = double(HRIR_L2702.irChTwo);
hrirGridAziRad = double(HRIR_L2702.azimuth.');
hrirGridZenRad = double(HRIR_L2702.elevation.'); % the elevation angles actually contain zenith data between 0..pi
fs = double(HRIR_L2702.fs);
clear HRIR_L2702;

[sig, sigFs] = audioread(sigPath);
assert(sigFs == fs, 'Mismatch in sampling frequencies.');

smaSig = fftfilt(srirSmaStruct.rir, sig);
emaSig = fftfilt(srirEmaStruct.rir, sig);

%% get filters for the LS, MagLS, eMagLS and eMagLS2 renderers
fprintf('Computing rendering filters \n');
% LS solution
[wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, order, shDefinition);

% MagLS
[wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    order, fs, filterLen, shDefinition);

% EMagLS
[wEMlsL, wEMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    smaRadius, smaGridAziRad, smaGridZenRad, order, fs, filterLen, shDefinition);

% EMagLS2
[wEMls2L, wEMls2R] = getEMagLs2Filters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    smaRadius, smaGridAziRad, smaGridZenRad, order, fs, filterLen, shDefinition);

% EMagLS for EMA using CHs
[wEMlsLEmaCh, wEMlsREmaCh] = getEMagLsFiltersEMAinCH(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    emaRadius, emaGridAziRad, order, fs, filterLen, shDefinition);

% EMagLS for EMA using SHs
[wEMlsLEmaSh, wEMlsREmaSh] = getEMagLsFiltersEMAinSH(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    emaRadius, emaGridAziRad, order, fs, filterLen, shDefinition);

%% SH transform and radial filter (only for LS and conventional MagLS)
EncSma = getSH(order, [smaGridAziRad, smaGridZenRad], shDefinition); % SMA encoder
EncEma = getCH(order, emaGridAziRad, shDefinition); % EMA encoder

shSigSma = smaSig * pinv(EncSma.');
chSigEma = emaSig * pinv(EncEma.');

J = getChToShExpansionMatrix(order, shDefinition);
shSigEma = chSigEma * J.';

% parameters for the radial filter
params.order = order;
params.fs = fs;
params.irLen = filterLen;
params.oversamplingFactor = 1;
params.radialFilter = 'tikhonov';
params.regulConst = 1e-2;
params.smaRadius = smaRadius;
params.smaDesignAziZenRad = [smaGridAziRad, smaGridZenRad];
params.waveModel = 'planeWave';
params.arrayType = 'rigid';
params.nfft = params.oversamplingFactor * params.irLen;

shSigSmaRadFiltered = applyRadialFilter(shSigSma, params);

%% render to binaural
fprintf('Rendering \n');

% the LS and MagLS renderers need the radial filtered signals as input
binLs = binauralDecode(shSigSmaRadFiltered, fs, wLsL, wLsR, fs);
binMls = binauralDecode(shSigSmaRadFiltered, fs, wMlsL, wMlsR, fs);

% the eMagLS renderer needs the unfiltered SH-domain signal as input
binEMls = binauralDecode(shSigSma, fs, wEMlsL, wEMlsR, fs);

% the eMagLS2 renderer needs the raw microphone signals as input
binEMls2 = binauralDecode(smaSig, fs, wEMls2L, wEMls2R, fs);

% for EMA rendering, there two options: using SHs (allowing for 3DoF
% headtracking) or using CHs (supporting 1DoF headtracking)
binEMlsEmaCh = binauralDecode(chSigEma, fs, wEMlsLEmaCh, wEMlsREmaCh, fs);
binEMlsEmaSh = binauralDecode(shSigEma, fs, wEMlsLEmaSh, wEMlsREmaSh, fs);

fprintf('Normalizing binaural renderings \n');
binLs = binLs ./ max(abs(binLs(:))) * 0.5;
binMls = binMls ./ max(abs(binMls(:))) * 0.5;
binEMls = binEMls ./ max(abs(binEMls(:))) * 0.5;
binEMls2 = binEMls2 ./ max(abs(binEMls2(:))) * 0.5;
binEMlsEmaCh = binEMlsEmaCh ./ max(abs(binEMlsEmaCh(:))) * 0.5;
binEMlsEmaSh = binEMlsEmaSh ./ max(abs(binEMlsEmaSh(:))) * 0.5;

fprintf('Starting playback, source should appear at %s° azimuth and %s° elevation ... \n', ...
    num2str(srirSmaStruct.dirSoundDoaAziEleRad(1)*180/pi), num2str(srirSmaStruct.dirSoundDoaAziEleRad(2)*180/pi));
fprintf('Playing back LS binaural rendering ... \n');
playblocking(audioplayer(binLs, fs));

pause(0.5);
fprintf('Playing back MagLS binaural rendering ... \n');
playblocking(audioplayer(binMls, fs));

pause(0.5);
fprintf('Playing back eMagLS binaural rendering ... \n');
playblocking(audioplayer(binEMls, fs));

pause(0.5);
fprintf('Playing back eMagLS2 binaural rendering ... \n');
playblocking(audioplayer(binEMls2, fs));

pause(0.5);
fprintf('Playing back eMagLS binaural rendering from EMA using CHs ... \n');
playblocking(audioplayer(binEMlsEmaCh, fs));

pause(0.5);
fprintf('Playing back eMagLS binaural rendering from EMA using SHs ... \n');
playblocking(audioplayer(binEMlsEmaSh, fs));