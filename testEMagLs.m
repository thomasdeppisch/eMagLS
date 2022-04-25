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

clear; clc; close all;

addpath(genpath('dependencies/'));

%% configuration
filterLen = 512;

[hrirFile, hrirUrl] = deal('resources/HRIR_L2702.mat', ...
    'https://zenodo.org/record/3928297/files/HRIR_L2702.mat');

[smaRecordingLength, smaRecordingFile, smaRecordingUrl] = deal( ...
    13.5, ... % in seconds
    'resources/Acappella_Eigenmike_Raw_32ch.wav', ...
    'https://zenodo.org/record/3477602/files/09%203D-MARCo%20Samples_Acappella.zip');
% define SMA geometry for recording file (here Eigenmike EM32)
shOrder       = 4;
micRadius     = 0.042; % in m
micGridAziRad = deg2rad([0;32;0;328;0;45;69;45;0;315;291;315;91;90;90;89;180;212;180;148;180;225;249;225;180;135;111;135;269;270;270;271]);
micGridZenRad = deg2rad([69;90;111;90;32;55;90;125;148;125;90;55;21;58;121;159;69;90;111;90;32;55;90;125;148;125;90;55;21;58;122;159]);

%% load data
tic; % start measuring execution time

fprintf('Downloading HRIR dataset ... ');
if isfile(hrirFile)
    fprintf('already exists ... skipped.\n');
else
    downloadAndExtractFile(hrirFile, hrirUrl);
end

fprintf('Loading file "%s" ... ', hrirFile);
load(hrirFile);
hL = double(HRIR_L2702.irChOne);
hR = double(HRIR_L2702.irChTwo);
hrirGridAziRad = double(HRIR_L2702.azimuth.');
hrirGridZenRad = double(HRIR_L2702.elevation.'); % the elevation angles actually contain zenith data between 0..pi
fs = double(HRIR_L2702.fs);
clear HRIR_L2702;
fprintf('done.\n');

fprintf('Downloading SMA recording ... ');
if isfile(smaRecordingFile)
    fprintf('already exists ... skipped.\n');
else
    downloadAndExtractFile(smaRecordingFile, smaRecordingUrl);
end

fprintf('Loading file "%s" ... ', smaRecordingFile);
[smaRecording, smaFs] = audioread(smaRecordingFile);
if smaFs ~= fs; error('Mismatch in sampling frequencies.'); end
if size(smaRecording, 1) / fs > smaRecordingLength
    fprintf('truncating length to %.1f s ... ', smaRecordingLength);
    smaRecording = smaRecording(1:fs * smaRecordingLength, :); % truncate
end
clear smaFs smaRecordingLength;
fprintf('done.\n\n');

%% get filters for the LS, MagLS, eMagLS and eMagLS2 renderers
fprintf('Computing LS rendering filters ... with %d samples ... ', filterLen);
[wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, shOrder);
fprintf('done.\n');

fprintf('Computing MagLS rendering filters ... with %d samples ... ', filterLen);
[wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    shOrder, fs, filterLen, false);
fprintf('done.\n');

fprintf('Computing eMagLS rendering filters ... with %d samples ... ', filterLen);
[wEMlsL, wEMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, micGridZenRad, shOrder, fs, filterLen, false);
fprintf('done.\n');

fprintf('Computing eMagLS2 rendering filters ... with %d samples ... ', filterLen);
[wEMls2L, wEMls2R] = getEMagLs2Filters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, micGridZenRad, fs, filterLen, false);
fprintf('done.\n\n');

%% SH transform and radial filter (for LS and conventional MagLS)
fprintf('Transforming recording into SH domain at N=%d ... ', shOrder);
E = getSH(shOrder, [micGridAziRad, micGridZenRad], 'real');
shRecording = smaRecording * pinv(E)';
fprintf('done.\n');

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

fprintf('Computing and applying radial filters ... ');
shRecordingRadFiltered = applyRadialFilter(shRecording, params);
fprintf('done.\n\n');

%% render binaurally
% the LS and MagLS renderers need the radial filtered signals as input
fprintf('Computing LS binaural rendering ... ');
lsBin = binauralDecode(shRecordingRadFiltered, fs, wLsL, wLsR, fs);
fprintf('done.\n');

fprintf('Computing MagLS binaural rendering ... ');
magLsBin = binauralDecode(shRecordingRadFiltered, fs, wMlsL, wMlsR, fs);
fprintf('done.\n');

fprintf('Computing eMagLS binaural rendering ... ');
% the eMagLS renderer needs the unfiltered SH-domain signal as input
eMagLsBin = binauralDecode(shRecording, fs, wEMlsL, wEMlsR, fs);
fprintf('done.\n');

fprintf('Computing eMagLS2 binaural rendering ... ');
% the eMagLS2 renderer needs the raw microphone signals as input
eMagLs2Bin = binauralDecode(smaRecording, fs, wEMls2L, wEMls2R, fs);
fprintf('done.\n');

fprintf('Normalizing binaural renderings ... ');
lsBin = lsBin ./ max(abs(lsBin(:))) * 0.5;
magLsBin = magLsBin ./ max(abs(magLsBin(:))) * 0.5;
eMagLsBin = eMagLsBin ./ max(abs(eMagLsBin(:))) * 0.5;
eMagLs2Bin = eMagLs2Bin ./ max(abs(eMagLs2Bin(:))) * 0.5;
fprintf('done.\n\n');

%% listen to binaural renderings
% especially concentrate on the very high frequencies to spot differences
% between MagLS and eMagLS (good headphones needed)
fprintf('Playing back LS rendering ... ');
playblocking(audioplayer(lsBin, fs));
fprintf('done.\n');

pause(0.5);
fprintf('Playing back MagLS rendering ... ');
playblocking(audioplayer(magLsBin, fs));
fprintf('done.\n');

pause(0.5);
fprintf('Playing back eMagLS rendering ... ');
playblocking(audioplayer(eMagLsBin, fs));
fprintf('done.\n');

pause(0.5);
fprintf('Playing back eMagLS2 rendering ... ');
playblocking(audioplayer(eMagLs2Bin, fs));
fprintf('done.\n\n');

%%
fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));
