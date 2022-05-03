function [smairMat, params] = getSMAIRMatrix(params)
% [smairMat, params] = getSMAIRMatrix(params)
% 
% calculates SMAIR transform matrices, i.e. SMA processing (plane-wave radiation,
% scattering, mic encoding, radial filtering) but without any sources
%
% smairMat      .. numShsSimulation x numShsOut x numFreqs
%               .. or (if returnRawMicSigs) numShsSimulation x numMics x numFreqs
% params        .. updated parameters
%
% Some parameters are limited to the following options:
% arrayType                     .. {'rigid', 'open'}
% radialFilter                  .. {'none', 'full', 'regul', 'softLimit', 'em32-zStyle', 'em32-zStyle-ffEq'}
% waveModel                     .. {'planeWave', 'pointSource'}
% simulateAliasing              .. {true, false}
% zStyleMaxRe                   .. {0, 1}
% replaceScatteringByRadFilt    .. {true, false} -> use radial filter as regularized scattering
%                                  effect to limit gain in equalization filter
% returnRawDiaphSigs            .. only return mic signals without SH transform at the output
% 
% Further parameters:
% smaDesignAziZenRad, order, simulationOrder, fs, smaRadius, sourceDist,
% noiseGainDb, oversamplingFactor, irLen, shDefinition, shFunction
%
% Most parameters have default options! (default is a plane-wave em32 simulation)
%
% This software is licensed under a Non-Commercial Software License 
% (see https://github.com/thomasdeppisch/eMagLS/blob/master/LICENSE for full details).
%
% Thomas Deppisch, 2021
    
    % parse params
    if (nargin < 1 || ~isfield(params,'smaDesignAziZenRad'))
        smaDesign = load('des.3.32.7.txt'); % t-design on which em32 is based (not perfectly the same layout)
        smaDesign = reshape(smaDesign,3,length(smaDesign)/3)';
        [micsAzi, micsEle, ~] = cart2sph(smaDesign(:,1),smaDesign(:,2),smaDesign(:,3));
        micsZen = pi/2 - micsEle;
        params.smaDesignAziZenRad = [micsAzi, micsZen];
    end
    if (nargin < 1 || ~isfield(params,'order'))
        params.order = 4;
    end
    if (nargin < 1 || ~isfield(params,'fs'))
        params.fs = 48000;
    end
    if (nargin < 1 || ~isfield(params,'smaRadius'))
        params.smaRadius = 0.042;
    end
    if (nargin < 1 || ~isfield(params,'arrayType'))
        params.arrayType = 'rigid';
    end
    if (nargin < 1 || ~isfield(params,'radialFilter'))
        params.radialFilter = 'regul';
    end
    if (nargin < 1 || ~isfield(params,'simulationOrder'))
        params.simulationOrder = 32;
    end
    if (nargin < 1 || ~isfield(params,'sourceDist'))
        params.sourceDist = 2;
    end
    if (nargin < 1 || ~isfield(params,'waveModel'))
        params.waveModel = 'planeWave';
    end
    if (nargin < 1 || ~isfield(params,'noiseGainDb'))
        params.noiseGainDb = 20;
    end
    if (nargin < 1 || ~isfield(params,'zStyleMaxRe'))
        params.zStyleMaxRe = 1;
    end
    if (nargin < 1 || ~isfield(params,'sourcePosCart'))
        params.sourcePosCart = [params.sourceDist;0;0];
    end
    if (nargin < 1 || ~isfield(params,'oversamplingFactor'))
        params.oversamplingFactor = 4;
    end
    if (nargin < 1 || ~isfield(params,'irLen'))
        params.irLen = 2048;
    end
    if (nargin < 1 || ~isfield(params,'returnRawMicSigs'))
        params.returnRawMicSigs = false;
    end
    if (nargin < 1 || ~isfield(params,'shDefinition'))
        params.shDefinition = 'real';
    end
    if (nargin < 1 || ~isfield(params,'shFunction'))
        params.shFunction = @getSH;
    end
    
    C = 343; % speed of sound in m/s

    nfft = params.oversamplingFactor*params.irLen;
    f = linspace(0,params.fs/2,nfft/2+1)';
    k = 2*pi*f/C;
    kr = k * params.smaRadius;
    params.sourceDist = norm(params.sourcePosCart); % if params.sourcePosCart is set this will overwrite the sourceDist setting!
    krSource = k * params.sourceDist;
    numShsOut = (params.order+1)^2;
    numShsSimulation = (params.simulationOrder+1)^2;
    numFreqs = length(k);
    numMics = size(params.smaDesignAziZenRad, 1);

    % set radial filtering for waveModel + arrayType combination
    switch lower(params.arrayType)
        case 'rigid'
            bn = @(N_,kr_) 1i ./ ((kr_.').^2 .* sph_besselh_diff(N_, kr_).');

        case 'open'
            bn = @(N_,kr_) sph_besselj(N_, kr_).';
            
        case 'directional' % open array with first-order directional mics -> see Politis, Array Response Simulator
            % 0 .. omni, 0.5 .. cardioid, 1 .. fig-of-eight
            bn = @(N_,kr_) (params.dirCoeff*sph_besselj(N_, kr_).' - 1i*(1-params.dirCoeff)*sph_besselj_diff(N_, kr_).');

        otherwise
            error('Unkown arrayType parameter "%s".', params.arrayType);
    end
    
    % include actual microphone processing to simulate aliasing
    n = (0:params.simulationOrder);
%     fprintf('with @%s("%s") ... ', func2str(params.shFunction), params.shDefinition);
    YMics = params.shFunction(params.simulationOrder, params.smaDesignAziZenRad, params.shDefinition).';
    YMicsLow = YMics(1:(params.order+1)^2, :);

    switch lower(params.waveModel)
        case 'planewave'
            bnAll = 4*pi*1i.^n .* bn(params.simulationOrder, kr).';

        case 'pointsource'
            hnAll = sph_besselh(params.simulationOrder, krSource);
            bnAll = 4*pi*(-1i) .* k .* hnAll .* bn(params.simulationOrder, kr).';

        otherwise
            error('Unkown waveModel parameter "%s".', params.waveModel);
    end

    pMics = zeros(numMics, numShsSimulation, numFreqs);

    for ii = 1:numFreqs
        Bn = diag(sh_repToOrder(bnAll(ii,:).').');
        pMics(:,:,ii) = YMics' * Bn;
    end

    pMics(isnan(pMics)) = 0; % remove singularity for f = 0

    pN = zeros(numShsOut, numShsSimulation, numFreqs);
    pinvYMicsLow = pinv(YMicsLow');
    for ii = 1:numFreqs
        pN(:,:,ii) = pinvYMicsLow * pMics(:,:,ii);
    end

    % apply radial filtering
    if (~strcmpi(params.radialFilter, 'none'))
        radFilts = getRadialFilter(params);

        smairMat = zeros(numShsOut, numShsSimulation, numFreqs);
        for ii = 1:numFreqs
            BnTi = diag(sh_repToOrder(radFilts(ii,:).').');
            smairMat(:,:,ii) = BnTi * pN(:,:,ii);
        end

        smairMat(isnan(smairMat)) = 0;
    else
        smairMat = pN;
    end

    if params.returnRawMicSigs
        smairMat = pMics;
    end
end
