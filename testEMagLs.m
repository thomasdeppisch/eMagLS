% This script demonstrates how the different eMagLS renderers can be used
% to render a SMA recording binaurally. It also provides the opportunity to
% listen and compare different renderers.
%
% For more information, please refer to
%   T. Deppisch, H. Helmholz, and J. Ahrens,
%   “End-to-End Magnitude Least Squares Binaural Rendering of Spherical Microphone Array Signals,”
%   in 2021 Immersive and 3D Audio: from Architecture to Automotive (I3DA), 2021, pp. 1–7. doi: 10.1109/I3DA48870.2021.9610864.
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Hannes Helmholz, 2023

clear; clc; close all;

addpath(genpath('dependencies/'));

%% configuration
filterLen = 512; % in samples
shDefinition = 'complex';

[hrirFile, hrirUrl] = deal('resources/HRIR_L2702.mat', ...
    'https://zenodo.org/record/3928297/files/HRIR_L2702.mat');

smaRecordingFile = 'resources/Acappella_Eigenmike_Raw_32ch_short.wav';
% define SMA geometry for recording file (here Eigenmike EM32)
shOrder       = 4;
micRadius     = 0.042; % in m
micGridAziRad = deg2rad([0;32;0;328;0;45;69;45;0;315;291;315;91;90;90;89;180;212;180;148;180;225;249;225;180;135;111;135;269;270;270;271]);
micGridZenRad = deg2rad([69;90;111;90;32;55;90;125;148;125;90;55;21;58;121;159;69;90;111;90;32;55;90;125;148;125;90;55;21;58;122;159]);

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

fprintf('Loading file "%s" ... \n', smaRecordingFile);
[smaRecording, smaFs] = audioread(smaRecordingFile);
assert(smaFs == fs, 'Mismatch in sampling frequencies.');
if exist('smaRecordingLength', 'var') && size(smaRecording, 1) / fs > smaRecordingLength
    fprintf('truncating length to %.1f s ... \n', smaRecordingLength);
    smaRecording = smaRecording(1:fs * smaRecordingLength, :); % truncate
end

%% get filters for the LS, MagLS, eMagLS and eMagLS2 renderers
fprintf('Computing rendering filters \n');
[wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, shOrder, shDefinition);
[wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    shOrder, fs, filterLen, shDefinition);
[wEMlsL, wEMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, micGridZenRad, shOrder, fs, filterLen, shDefinition);

% EMA
[wEMlsLEma, wEMlsREma] = getEMagLsFiltersEma(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, shOrder, fs, filterLen, shDefinition);

% EMagLS2
[wEMls2L, wEMls2R] = getEMagLs2Filters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, micGridZenRad, shOrder, fs, filterLen, shDefinition);

%% SH transform and radial filter (for LS and conventional MagLS)
E = getSH(shOrder, [micGridAziRad, micGridZenRad], shDefinition).';
shRecording = smaRecording * pinv(E);

% parameters for the radial filter
params.order = shOrder;
params.fs = fs;
params.irLen = filterLen;
params.oversamplingFactor = 1;
params.radialFilter = 'tikhonov';
params.regulConst = 1e-2;
params.smaRadius = micRadius;
params.smaDesignAziZenRad = [micGridAziRad, micGridZenRad];
params.waveModel = 'planeWave';
params.arrayType = 'rigid';
params.nfft = params.oversamplingFactor * params.irLen;

shRecordingRadFiltered = applyRadialFilter(shRecording, params);

%% render binaurally
fprintf('Rendering \n');
% the LS and MagLS renderers need the radial filtered signals as input
binLs = binauralDecode(shRecordingRadFiltered, fs, wLsL, wLsR, fs);
binMls = binauralDecode(shRecordingRadFiltered, fs, wMlsL, wMlsR, fs);

% the eMagLS renderer needs the unfiltered SH-domain signal as input
binEMls = binauralDecode(shRecording, fs, wEMlsL, wEMlsR, fs);

% the eMagLS2 renderer needs the raw microphone signals as input
binEMls2 = binauralDecode(smaRecording, fs, wEMls2L, wEMls2R, fs);

fprintf('Normalizing binaural renderings \n');
binLs = binLs ./ max(abs(binLs(:))) * 0.5;
binMls = binMls ./ max(abs(binMls(:))) * 0.5;
binEMls = binEMls ./ max(abs(binEMls(:))) * 0.5;
binEMls2 = binEMls2 ./ max(abs(binEMls2(:))) * 0.5;

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
