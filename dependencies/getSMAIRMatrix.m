function [smairMat, params] = getSMAIRMatrix(params)
    % returns the SMAIR transform matrices, i.e. SMA processing (plane-wave
    % radiation, scattering, mic encoding, radial filtering) but without
    % any sources
    %
    % smairMat .. numShsSimulation x numShsOut x numFreqs
    %          .. or (if returnRawMicSigs) numShsSimulation x numMics x numFreqs
    %
    % some params are limited to the following options:
    % arrayType: {'rigid', 'open'}
    % radialFilter: {'none', 'full', 'regul', 'softLimit', 'em32-zStyle', 'em32-zStyle-ffEq'}
    % waveModel: {'planeWave', 'pointSource'}
    % simulateAliasing: {true, false}
    % zStyleMaxRe: {0, 1}
    % replaceScatteringByRadFilt = {true, false} -> use radial filter as
    % regularized scattering effect to limit gain in equalization filter
    % returnRawDiaphSigs .. only return mic signals without SH transform at
    %                       the output
    % 
    % further params:
    % smaDesignAziZenRad, order, simulationOrder, fs, smaRadius,
    % sourceDist, noiseGainDb, oversamplingFactor, irLen
    %
    % most parameters have default options! (default is a plane-wave em32
    % simulation)
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
    if (nargin < 1 || ~isfield(params,'dirCoeff'))
        params.dirCoeff = 0;
    end
    
    c = 343;
    nfft = params.oversamplingFactor*params.irLen;
    f = linspace(0,params.fs/2,nfft/2+1)';
    k = 2*pi*f/c;
    kr = k * params.smaRadius;
    params.sourceDist = norm(params.sourcePosCart); % if params.sourcePosCart is set this will overwrite the sourceDist setting!
    %krSource = k * params.sourceDist;
    numShsOut = (params.order+1)^2;
    numShsSimulation = (params.simulationOrder+1)^2;
    numFreqs = length(k);
    numMics = size(params.smaDesignAziZenRad, 1);
    
    % include actual microphone processing to simulate aliasing
    YMics = getSH(params.simulationOrder, params.smaDesignAziZenRad, params.shDefinition).';    
    YMicsLow = YMics(1:(params.order+1)^2, :);

    bnAll = sphModalCoeffs(params.simulationOrder, kr.', params.arrayType, params.dirCoeff); % todo: add spherical waves
    bnAll(end,:) = real(bnAll(end,:));

    pMics = zeros(numMics, numShsSimulation, nfft);

    for ii = 1:numFreqs
        Bn = diag(sh_repToOrder(bnAll(ii,:).').');
        pMics(:,:,ii) = YMics.' * Bn;

        if ii > 1 && ii < numFreqs
            pMics(:,:,end-ii+2) = YMics.' * conj(Bn);
        end
    end

    pMics(isnan(pMics)) = 0; % remove singularity for f = 0

    pN = zeros(numShsOut, numShsSimulation, nfft);
    pinvYMicsLow = pinv(YMicsLow.');
    for ii = 1:nfft
        pN(:,:,ii) = pinvYMicsLow * pMics(:,:,ii);
    end

    % apply radial filtering
    if (~strcmp(params.radialFilter, 'none'))
        radFilts = getRadialFilter(params);

        smairMat = zeros(numShsOut, numShsSimulation, nfft);
        for ii = 1:numFreqs
            BnTi = diag(sh_repToOrder(radFilts(ii,:).').');
            smairMat(:,:,ii) = BnTi * pN(:,:,ii);

            if ii > 1 && ii < numFreqs
                smairMat(:,:,end-ii+2) = conj(BnTi) * pN(:,:,ii);
            end
        end

        smairMat(isnan(smairMat)) = 0;
    else
        smairMat = pN;
    end

    if params.returnRawMicSigs
        smairMat = pMics;
    end
    