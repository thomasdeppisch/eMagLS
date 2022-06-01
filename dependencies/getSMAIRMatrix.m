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
% noiseGainDb, oversamplingFactor, irLen, dirCoeff, shDefinition, shFunction
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
    if (nargin < 1 || ~isfield(params,'dirCoeff'))
        params.dirCoeff = 0;
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

    nfft = params.oversamplingFactor * params.irLen;
    f = linspace(0, params.fs/2, nfft/2+1).';
    params.sourceDist = norm(params.sourcePosCart); % if params.sourcePosCart is set this will overwrite the sourceDist setting!
    numShsOut = (params.order+1)^2;
    numShsSimulation = (params.simulationOrder+1)^2;
    numPosFreqs = length(f);
    numMics = size(params.smaDesignAziZenRad, 1);

    % include actual microphone processing to simulate aliasing
    Y_Hi = params.shFunction(params.simulationOrder, params.smaDesignAziZenRad, params.shDefinition);
    Y_Lo_pinv = pinv(Y_Hi(:, 1:numShsOut));

    % set radial filtering for waveModel + arrayType combination
    if strcmpi(params.waveModel, 'pointSource')
        % TODO: Add spherical waves
        error('WaveModel parameter "%s" not yet implemented.', params.waveModel);
    end
    % TODO: It is not quite clear why there is a minus required here for
    % the rendered BRIRs to start with a positive peak (which seems
    % reasonable). Without the minus, the resulting BRIRs are inverted.
    bnAll = -sphModalCoeffs(params.simulationOrder, 2*pi*f/C * params.smaRadius, ...
        params.arrayType, params.dirCoeff).';

    pMics = zeros(numMics, numShsSimulation, nfft);
    pN = zeros(numShsOut, numShsSimulation, nfft);
    for k = 1:numPosFreqs
        Bn = diag(sh_repToOrder(bnAll(:,k)));
        pMics(:,:,k) = Y_Hi * Bn;
        if k == numPosFreqs && ~mod(nfft, 2) % is even
            pMics(:,:,k) = Y_Hi * real(Bn);
        elseif ~isreal(Y_Hi) && k > 1
            k_neg = nfft-k+2;
            pMics(:,:,k_neg) = Y_Hi * conj(Bn);
        end

        if ~params.returnRawMicSigs
            pN(:,:,k) = Y_Lo_pinv * pMics(:,:,k);
            if k == numPosFreqs && ~mod(nfft, 2) % is even
                pN(:,:,k) = Y_Lo_pinv * real(pMics(:,:,k));
            elseif ~isreal(Y_Hi) && k > 1
                pN(:,:,k_neg) = Y_Lo_pinv * pMics(:,:,k_neg);
            end
        end
    end

    if params.returnRawMicSigs
        smairMat = pMics;
    else
        smairMat = pN;
        
        if ~strcmpi(params.radialFilter, 'none')
            % apply radial filtering
            radFilts = getRadialFilter(params).';
            for k = 1:numPosFreqs
                BnTi = diag(sh_repToOrder(radFilts(:,k)));
                smairMat(:,:,k) = BnTi * smairMat(:,:,k);
                if k == numPosFreqs && ~mod(nfft, 2) % is even
                    smairMat(:,:,k) = real(BnTi) * smairMat(:,:,k);
                elseif ~isreal(Y_Hi) && k > 1
                    k_neg = nfft-k+2;
                    smairMat(:,:,k_neg) = conj(BnTi) * smairMat(:,:,k_neg);
                end
            end
        end
    end
end
