% This script demonstrates how the different eMagLS renderers can be used
% to render a SMA recording binaurally.
% 
% Thomas Deppisch, 2021

clear all
close all

addpath(genpath('dependencies/'))
[smaRecording, ~] = audioread('resources/EM32ch_hall_social1_short.wav');

%%
shOrder = 4;
filterLen = 512;

% define SMA design (em32)
micRadius = 0.042;
micGridAziRad = pi/180 * [0;32;0;328;0;45;69;45;0;315;291;315;91;90;90;89;180;212;180;148;180;225;249;225;180;135;111;135;269;270;270;271];
micGridZenRad = pi/180 * [69;90;111;90;32;55;90;125;148;125;90;55;21;58;121;159;69;90;111;90;32;55;90;125;148;125;90;55;21;58;122;159];

% load HRIR set
load HRIR_L2702.mat
hL = double(HRIR_L2702.irChOne);
hR = double(HRIR_L2702.irChTwo);
hrirGridAziRad = HRIR_L2702.azimuth';
hrirGridZenRad = HRIR_L2702.elevation'; % the elevation angles actually contain zenith data between 0..pi
fs = double(HRIR_L2702.fs);

[wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, shOrder);
[wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, shOrder, fs, filterLen, false);
[wEMlsL, wEMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, micRadius, micGridAziRad, micGridZenRad, shOrder, fs, filterLen, false);
[wEMls2L, wEMls2R] = getEMagLs2Filters(hL, hR, hrirGridAziRad, hrirGridZenRad, micRadius, micGridAziRad, micGridZenRad, fs, filterLen, false);

%% sh transform and radial filter (for LS and conventional MagLS)
E = getSH(shOrder, [micGridAziRad, micGridZenRad], 'real');
shRecording = smaRecording * pinv(E)';

params.order = shOrder;
params.fs = fs;
params.irLen = 256;
params.oversamplingFactor = 1;
params.radialFilter = 'regul';
params.smaRadius = micRadius;
params.smaDesignAziZenRad = [micGridAziRad, micGridZenRad];
params.waveModel = 'planeWave';
params.arrayType = 'rigid';
params.nfft = params.oversamplingFactor*params.irLen;

shRecordingRadFiltered = applyRadialFilter(shRecording, params);

%% render binaurally
% the LS and MagLS renderers need the radial filtered signals as input
lsBin = binauralDecode(shRecordingRadFiltered, fs, wLsL, wLsR, fs);
magLsBin = binauralDecode(shRecordingRadFiltered, fs, wMlsL, wMlsR, fs);

% the eMagLS renderer needs the unfiltered SH-domain signal
eMagLsBin = binauralDecode(shRecording, fs, wEMlsL, wEMlsR, fs);

% the eMagLS2 renderer needs the raw microphone signals
eMagLs2Bin = binauralDecode(smaRecording, fs, wEMls2L, wEMls2R, fs);

%% listen to binaural renderings
% LS
sound(lsBin./max(abs(lsBin(:))) * 0.5, fs);
pause(7)
% MagLS
sound(magLsBin./max(abs(magLsBin(:))) * 0.5, fs);
pause(7)
% eMagLS
sound(eMagLsBin./max(abs(eMagLsBin(:))) * 0.5, fs);
pause(7)
% eMagLS2
sound(eMagLs2Bin./max(abs(eMagLs2Bin(:))) * 0.5, fs);