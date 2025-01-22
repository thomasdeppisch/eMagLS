% This script demonstrates how eMagLS can be used with arbitrary arrays
% based on measured array transfer functions (ATFs).
% This example renders signals from a pair of smartglasses.
%
% For more information, please refer to
%   Thomas Deppisch, Nils Meyer-Kahlen, Sebastia Amengual Gari, 
%   "Blind Identification of Binaural Room Impulse Responses from Smart Glasses", 
%   IEEE/ACM Transactions on Audio, Speech, and Language Processing, vol. 32, pp. 4052-4065, 2024, 
%   doi: 10.1109/TASLP.2024.3454964.
%
% This software is licensed under a Non-Commercial Software License
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Thomas Deppisch, 2025

clear; clc; close all;

addpath(genpath('dependencies/'));
addpath(genpath('lib/'));

%% load
% load HRIRs
[hrirFile, hrirUrl] = deal('./resources/HRIR_L2702.mat', ...
    'https://zenodo.org/record/3928297/files/HRIR_L2702.mat');

% load a speech signal
[sig,sigFs] = audioread('./resources/decemberTour.wav'); 

% load a RIR
rirStruct = load('./resources/meetingRoom_leftLsp.mat');
rirFs = rirStruct.fs;

% load ATFs
atfStructFullSphere = load('./resources/glasses_on_HATS_ATFs_sphere.mat');
atfIrGridFullSphere = atfStructFullSphere.atfIrs;
atfAziGridRadFS = atfStructFullSphere.atfGridAziEleDeg(:,1) * pi/180;
atfZenGridRadFS = (90-atfStructFullSphere.atfGridAziEleDeg(:,2)) * pi/180;
atfFs = atfStructFullSphere.fs;

fprintf('Downloading HRIR dataset ... \n');
if isfile(hrirFile)
    fprintf('already exists ... skipped.\n');
else
    downloadAndExtractFile(hrirFile, hrirUrl);
end

% load HRIR
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
hrirFs = double(HRIR_L2702.fs);
clear HRIR_L2702;

assert(hrirFs == sigFs && hrirFs == atfFs && hrirFs == rirFs, 'Sample Rate Mismatch!')

%% create array signal
arraySig = fftfilt(rirStruct.roomIRs, sig);

%% calculate eMagLS rendering filters
eMagLsFilterLen = 256;
eMagLsCutOnFreq = 2000;
[wMlsL, wMlsR] = getEMagLsFiltersFromAtf(hL, hR, [hrirGridAziRad, hrirGridZenRad], atfIrGridFullSphere, ...
                                        [atfAziGridRadFS atfZenGridRadFS], sigFs, eMagLsFilterLen, eMagLsCutOnFreq);

%% render to binaural
fprintf('Rendering \n');

binEMls = binauralDecode(arraySig, hrirFs, wMlsL, wMlsR, sigFs);

% normalize
binEMls = binEMls ./ max(abs(binEMls(:))) * 0.5;

fprintf('Playing back eMagLS binaural rendering ... \n');
playblocking(audioplayer(binEMls, sigFs));

