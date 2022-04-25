% This script demonstrates how the different eMagLS renderers can be used
% to render a SMA recording binaurally. It also provides the opportunity to
% listen and compare different renderers.
%
% For more information about the renderer, please refer to 
% T. Deppisch, H. Helmholz, J. Ahrens, "End-to-End Magnitude Least Squares Binaural Rendering 
% of Spherical Microphone Array Signals," International 3D Audio Conference (I3DA), 2021.
% 
% This software is licensed under a Non-Commercial Software License 
% (see https://github.com/thomasdeppisch/eMagLS/blob/master/LICENSE for full details).
%
% Thomas Deppisch, 2021

clear all
close all

addpath(genpath('dependencies/'));

%% get filters for the LS, MagLS, eMagLS and eMagLS2 renderers
shOrder = 4;
filterLen = 512;

[hrirFile, hrirUrl] = deal('resources/HRIR_L2702.mat', ...
    'https://zenodo.org/record/3928297/files/HRIR_L2702.mat');

[smaRecordingFile, smaRecordingUrl] = deal('resources/Acappella_Eigenmike_Raw_32ch.wav', ...
    'https://zenodo.org/record/3477602/files/09%203D-MARCo%20Samples_Acappella.zip');

% define SMA design (em32)
micRadius = 0.042;
micGridAziRad = pi/180 * [0;32;0;328;0;45;69;45;0;315;291;315;91;90;90;89;180;212;180;148;180;225;249;225;180;135;111;135;269;270;270;271];
micGridZenRad = pi/180 * [69;90;111;90;32;55;90;125;148;125;90;55;21;58;121;159;69;90;111;90;32;55;90;125;148;125;90;55;21;58;122;159];

%% get filters for the LS, MagLS, eMagLS and eMagLS2 renderers
fprintf('Downloading HRIR dataset ... ');
if isfile(hrirFile)
    fprintf('already exists ... skipped.\n');
else
    downloadAndExtractFile(hrirFile, hrirUrl);
end

% load HRIR set
load(hrirFile);
hL = double(HRIR_L2702.irChOne);
hR = double(HRIR_L2702.irChTwo);
hrirGridAziRad = double(HRIR_L2702.azimuth.');
hrirGridZenRad = double(HRIR_L2702.elevation.'); % the elevation angles actually contain zenith data between 0..pi
fs = double(HRIR_L2702.fs);

% retrieve rendering filters
[wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, shOrder);
[wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, shOrder, fs, filterLen, false);
[wEMlsL, wEMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, micRadius, micGridAziRad, micGridZenRad, shOrder, fs, filterLen, false);
[wEMls2L, wEMls2R] = getEMagLs2Filters(hL, hR, hrirGridAziRad, hrirGridZenRad, micRadius, micGridAziRad, micGridZenRad, fs, filterLen, false);

%% SH transform and radial filter (for LS and conventional MagLS)
fprintf('Downloading SMA recording ... ');
if isfile(smaRecordingFile)
    fprintf('already exists ... skipped.\n');
else
    downloadAndExtractFile(smaRecordingFile, smaRecordingUrl);
end

% load SMA recording
[smaRecording, smaFs] = audioread(smaRecordingFile);
if smaFs ~= fs; error('Mismatch in sampling frequencies.'); end

E = getSH(shOrder, [micGridAziRad, micGridZenRad], 'real');
shRecording = smaRecording * pinv(E)';

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
params.nfft = params.oversamplingFactor*params.irLen;

% apply radial filter
shRecordingRadFiltered = applyRadialFilter(shRecording, params);

%% render binaurally
% the LS and MagLS renderers need the radial filtered signals as input
lsBin = binauralDecode(shRecordingRadFiltered, fs, wLsL, wLsR, fs);
magLsBin = binauralDecode(shRecordingRadFiltered, fs, wMlsL, wMlsR, fs);

% the eMagLS renderer needs the unfiltered SH-domain signal as input
eMagLsBin = binauralDecode(shRecording, fs, wEMlsL, wEMlsR, fs);

% the eMagLS2 renderer needs the raw microphone signals as input
eMagLs2Bin = binauralDecode(smaRecording, fs, wEMls2L, wEMls2R, fs);

%% listen to binaural renderings
% especially concentrate on the very high frequencies to spot differences
% between MagLS and eMagLS (good headphones needed)

% LS
disp('Playing back LS rendering')
sound(lsBin./max(abs(lsBin(:))) * 0.5, fs);
pause(8)
% MagLS
disp('Playing back MagLS rendering')
sound(magLsBin./max(abs(magLsBin(:))) * 0.5, fs);
pause(8)
% eMagLS
disp('Playing back eMagLS rendering')
sound(eMagLsBin./max(abs(eMagLsBin(:))) * 0.5, fs);
pause(8)
% eMagLS2
disp('Playing back eMagLS2 rendering')
sound(eMagLs2Bin./max(abs(eMagLs2Bin(:))) * 0.5, fs);
