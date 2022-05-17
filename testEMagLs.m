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
applyDiffusenessConst = false; % for MagLS, eMagLS and eMagLS2
shDefinition = 'real'; % or e.g. 'complex' # TODO: Fix complex rendering

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

global DO_VERIFY_REFERENCE DO_OVERRIDE_REFERENCE DO_EXPORT_RENDERING DO_PLAYBACK_RENDERING
DO_VERIFY_REFERENCE   = true;
DO_OVERRIDE_REFERENCE = false;
DO_EXPORT_RENDERING   = true;
DO_PLAYBACK_RENDERING = true;

% DO_VERIFY_REFERENCE   = false;
% DO_OVERRIDE_REFERENCE = true;
% DO_EXPORT_RENDERING   = false;
% DO_PLAY_RENDERING     = false;

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
assert(smaFs == fs, 'Mismatch in sampling frequencies.');
if size(smaRecording, 1) / fs > smaRecordingLength
    fprintf('truncating length to %.1f s ... ', smaRecordingLength);
    smaRecording = smaRecording(1:fs * smaRecordingLength, :); % truncate
end
clear smaFs smaRecordingLength;
fprintf('done.\n\n');

%% get filters for the LS, MagLS, eMagLS and eMagLS2 renderers
fprintf('Computing LS rendering filters ... with %d samples ... ', size(hL, 1));
[wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, shOrder, shDefinition);
fprintf('done.\n');

fprintf('Computing MagLS rendering filters ... with %d samples ... ', filterLen);
[wMlsL, wMlsR] = getMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    shOrder, fs, filterLen, applyDiffusenessConst, shDefinition);
fprintf('done.\n');

fprintf('Computing eMagLS rendering filters ... with %d samples ... ', filterLen);
[wEMlsL, wEMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, micGridZenRad, shOrder, fs, filterLen, ...
    applyDiffusenessConst, shDefinition);
fprintf('done.\n');

% % An alternative version which uses a different SH basis convention and implementation
% fprintf('Computing eMagLS rendering filters ... with %d samples ... ', filterLen);
% [wEMlsL, wEMlsR] = getEMagLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
%     micRadius, micGridAziRad, micGridZenRad, shOrder, fs, filterLen, ...
%     applyDiffusenessConst, shDefinition, @getSH_SFS);
% fprintf('done.\n');

fprintf('Computing eMagLS2 rendering filters ... with %d samples ... ', filterLen);
[wEMls2L, wEMls2R] = getEMagLs2Filters(hL, hR, hrirGridAziRad, hrirGridZenRad, ...
    micRadius, micGridAziRad, micGridZenRad, fs, filterLen, ...
    applyDiffusenessConst, shDefinition);
fprintf('done.\n\n');

%% verify rendering filters against provided reference
[hrirPath, refFiles, ~] = fileparts(hrirFile);
refStr = sprintf('%s_%dsamples_%dchannels_sh%d_%s_%%s', ...
    refFiles, filterLen, size(micGridAziRad, 1), shOrder, shDefinition);
refFiles = fullfile(hrirPath, [refStr, '.mat']);
clear hrirPath;

if DO_VERIFY_REFERENCE
    if DO_OVERRIDE_REFERENCE
        warning('Veryfying rendering filters ... skipped.');
    else
        % TODO: This verification could also check the match of other parameters

        refFile = sprintf(refFiles, 'LS');
        fprintf('Verifying LS rendering filters against "%s" ... ', refFile);
        ref = load(refFile);
        assertAllClose(wLsL, ref.wLsL);
        assertAllClose(wLsR, ref.wLsR);
        clear refFile ref;
        fprintf('done.\n');
    
        refFile = sprintf(refFiles, 'MagLS');
        fprintf('Verifying LS rendering filters against "%s" ... ', refFile);
        ref = load(refFile);
        assertAllClose(wMlsL, ref.wMlsL);
        assertAllClose(wMlsR, ref.wMlsR);
        clear refFile ref;
        fprintf('done.\n');
    
        refFile = sprintf(refFiles, 'eMagLS');
        fprintf('Verifying LS rendering filters against "%s" ... ', refFile);
        ref = load(refFile);
        assertAllClose(wEMlsL, ref.wEMlsL);
        assertAllClose(wEMlsR, ref.wEMlsR);
        clear refFile ref;
        fprintf('done.\n');
    
        refFile = sprintf(refFiles, 'eMagLS2');
        fprintf('Verifying LS rendering filters against "%s" ... ', refFile);
        ref = load(refFile);
        assertAllClose(wEMls2L, ref.wEMls2L);
        assertAllClose(wEMls2R, ref.wEMls2R);
        clear refFile ref;
        fprintf('done.\n\n');
    end
end

%% replace reference filters
if DO_OVERRIDE_REFERENCE
    refFile = sprintf(refFiles, 'LS'); %#ok<*UNRCH> 
    fprintf('Exporting LS rendering filters to "%s" ... ', refFile);
    save(refFile, 'wLsL', 'wLsR', 'hrirGridAziRad', 'hrirGridZenRad', 'shOrder', '-v7');
    fprintf('done.\n');
    
    refFile = sprintf(refFiles, 'MagLS');
    fprintf('Exporting MagLS rendering filters to "%s" ... ', refFile);
    save(refFile, 'wMlsL', 'wMlsR', 'hrirGridAziRad', 'hrirGridZenRad', ...
        'shOrder', 'fs', 'filterLen', 'applyDiffusenessConst', '-v7');
    fprintf('done.\n');
    
    refFile = sprintf(refFiles, 'eMagLS');
    fprintf('Exporting eMagLS rendering filters to "%s" ... ', refFile);
    save(refFile, 'wEMlsL', 'wEMlsR', 'hrirGridAziRad', 'hrirGridZenRad', ...
        'micRadius', 'micGridAziRad', 'micGridZenRad', ...
        'shOrder', 'fs', 'filterLen', 'applyDiffusenessConst', '-v7');
    fprintf('done.\n');
    
    refFile = sprintf(refFiles, 'eMagLS2');
    fprintf('Exporting eMagLS2 rendering filters to "%s" ... ', refFile);
    save(refFile, 'wEMls2L', 'wEMls2R', 'hrirGridAziRad', 'hrirGridZenRad', ...
        'micRadius', 'micGridAziRad', 'micGridZenRad', ...
        'fs', 'filterLen', 'applyDiffusenessConst', '-v7');
    fprintf('done.\n\n');
    clear refFile;
end

%% SH transform and radial filter (for LS and conventional MagLS)
fprintf('Transforming recording into SH domain at N=%d ... ', shOrder);
% This has to be adapted in case a different SH implementation is used
E = getSH(shOrder, [micGridAziRad, micGridZenRad], shDefinition);
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

%% export binaural renderings
[~, smaRecordingFile, ~] = fileparts(smaRecordingFile);
if DO_EXPORT_RENDERING
    audioFiles = sprintf('%s_%s.wav', smaRecordingFile, refStr);

    audioFile = sprintf(audioFiles, 'LS');
    fprintf('Exporting LS binaural rendering to "%s" ... ', audioFile);
    audiowrite(audioFile, lsBin, fs);
    fprintf('done.\n');

    audioFile = sprintf(audioFiles, 'MagLS');
    fprintf('Exporting LS binaural rendering to "%s" ... ', audioFile);
    audiowrite(audioFile, magLsBin, fs);
    fprintf('done.\n');

    audioFile = sprintf(audioFiles, 'eMagLS');
    fprintf('Exporting LS binaural rendering to "%s" ... ', audioFile);
    audiowrite(audioFile, eMagLsBin, fs);
    fprintf('done.\n');

    audioFile = sprintf(audioFiles, 'eMagLS2');
    fprintf('Exporting LS binaural rendering to "%s" ... ', audioFile);
    audiowrite(audioFile, eMagLs2Bin, fs);
    fprintf('done.\n\n');

    clear audioFiles audioFile;
end

%% listen to binaural renderings
% especially concentrate on the very high frequencies to spot differences
% between MagLS and eMagLS (good headphones needed)
if DO_PLAYBACK_RENDERING
    fprintf('Playing back LS binaural rendering ... ');
    playblocking(audioplayer(lsBin, fs));
    fprintf('done.\n');
    
    pause(0.5);
    fprintf('Playing back MagLS binaural rendering ... ');
    playblocking(audioplayer(magLsBin, fs));
    fprintf('done.\n');
    
    pause(0.5);
    fprintf('Playing back eMagLS binaural rendering ... ');
    playblocking(audioplayer(eMagLsBin, fs));
    fprintf('done.\n');
    
    pause(0.5);
    fprintf('Playing back eMagLS2 binaural rendering ... ');
    playblocking(audioplayer(eMagLs2Bin, fs));
    fprintf('done.\n\n');
end

%%
fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));

%% helper functions
% function Y = getSH_SHT(order, gridAziZenRad, shDefinition)
% % This is identical to the default implementation being used in the toolbox
%     % from Spherical-Harmonic-Transform toolbox
%     % $ git clone https://github.com/polarch/Spherical-Harmonic-Transform.git
%     Y = getSH(order, gridAziZenRad, shDefinition);
% end

% function Y = getSH_AKT(order, gridAziZenRad, shDefinition)
% % This uses a different SH implementation where this function has to match
% % the signature (parameters and output format) of `getSH()`
%     % from AKtools toolbox (run AKtoolsStart.m)
%     % $ svn checkout https://svn.ak.tu-berlin.de/svn/AKtools --username aktools --password ak
%     Y = AKsh(order, [], rad2deg(gridAziZenRad(:, 1)), ...
%         rad2deg(gridAziZenRad(:, 2)), shDefinition);
% end

% function Y = getSH_SFS(order, gridAziZenRad, shDefinition)
% % This uses a different SH implementation where this function has to match
% % the signature (parameters and output format) of `getSH()`
%     % from soundfieldsynthesis "Common" scripts
%     % $ git clone https://github.com/JensAhrens/soundfieldsynthesis.git
%     Y = zeros(size(gridAziZenRad, 1), (order+1)^2);
%     for n = 0 : order
%         for m = -n : n
%             Y(:, n^2+n+m+1) = sphharm(n, m, ...
%                 gridAziZenRad(:, 2), gridAziZenRad(:, 1), shDefinition);
%         end
%     end
% end

function assertAllClose(x1, x2, norm_tol)
    if nargin < 3; norm_tol = 1e-13; end

    norm_diff = max(abs(x1 - x2), [], 'all') / max(abs([x1, x2]), [], 'all');
    if norm_diff == 0
        fprintf('no difference ... ');
    elseif norm_diff < norm_tol
        fprintf('%.3g max. norm. abs. difference ... ', norm_diff);
    else
        error('Maximum normalized absolute mismatch of %.3g (greater than %.3g tolarance).', ...
            norm_diff, norm_tol);
    end
end
