function radFilts = getRadialFilter(params)
% return radial filters
% 
% Thomas Deppisch, 2021

    if (nargin < 1 || ~isfield(params,'radialFilter'))
        params.radialFilter = 'tikhonov';
    end
    if (nargin < 1 || ~isfield(params,'waveModel'))
        params.waveModel = 'planeWave';
    end
    if (nargin < 1 || ~isfield(params,'oversamplingFactor'))
        params.oversamplingFactor = 2;
    end
    if (nargin < 1 || ~isfield(params,'irLen'))
        params.irLen = 256;
    end
    if (nargin < 1 || ~isfield(params,'dirCoeff'))
            params.dirCoeff = 0;
    end

    C = 343; % speed of sound in m/s

    nfft = params.oversamplingFactor * params.irLen;
    f = linspace(0, params.fs/2, nfft/2+1).';

    if strcmpi(params.radialFilter, 'none')
        % nothing to do
        radFilts = ones(nfft/2+1, params.order+1);
        return;
    end

    % set radial filtering for waveModel + arrayType combination
    if strcmpi(params.waveModel, 'pointSource')
        % TODO: Add spherical waves
        error('WaveModel parameter "%s" not yet implemented.', params.waveModel);
    end
    bnAll = sphModalCoeffs(params.order, 2*pi*f/C * params.smaRadius, ...
        params.arrayType, params.dirCoeff);

    switch lower(params.radialFilter)
        case 'tikhonov'
            % tikhonov regularization,
            % e.g. used in Herzog, Habets, Eigenbeam-ESPRIT for DOA-Vector estimation
            if ~isfield(params,'regulConst')
                params.regulConst = 1e-2;
            end
            radFilts = conj(bnAll(:,1:params.order+1)) ...
                ./ (conj(bnAll(:,1:params.order+1)) .* bnAll(:,1:params.order+1) + params.regulConst);

        case 'softlimit'
            % Bernschuetz et al., Soft-Limiting der modalen Amplitudenverstärkung
            gainLimLin = 10^(params.noiseGainDb / 20);
            radFilts = 2*gainLimLin/pi * abs(bnAll(:,1:params.order+1)) ...
                ./ bnAll(:,1:params.order+1) ...
                .* atan(pi ./ (2*gainLimLin*abs(bnAll(:,1:params.order+1))));

        case 'full'
            radFilts = 1 ./ bnAll(:,1:params.order+1);

        otherwise
            error('Unkown radialFilter parameter "%s".', params.radialFilter);

    end

    if ~mod(nfft, 2) % is even
        radFilts(end, :) = abs(radFilts(end, :)); % Nyquist bin
    end
end
