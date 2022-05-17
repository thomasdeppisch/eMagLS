function radFilts = getRadialFilter(params)
% return radial filters
% Thomas Deppisch, 2021

if (nargin < 1 || ~isfield(params,'dirCoeff'))
        params.dirCoeff = 0;
end

c = 343;
nfft = params.oversamplingFactor*params.irLen;
f = linspace(0,params.fs/2,nfft/2+1)';
k = 2*pi*f/c;
kr = k * params.smaRadius;

% if strcmp(params.waveModel, 'pointSource')
%     krSource = k * params.sourceDist; % in this case params.sourceDist needs to be set by the user
% end

% set radial filtering for waveModel + arrayType combination
bnAll = sphModalCoeffs(params.order, kr.', params.arrayType, params.dirCoeff); % todo: add spherical waves
bnAll(end,:) = real(bnAll(end,:));
    
switch params.radialFilter
    case 'none'
        % nothing to do
        radFilts = ones(nfft/2+1, params.order+1);

    case 'tikhonov'
        % tikhonov regularization, 
        % e.g. used in Herzog, Habets, Eigenbeam-ESPRIT for DOA-Vector estimation
        if ~isfield(params,'regulConst')
            params.regulConst = 1e-2;
        end
        radFilts = (conj(bnAll(:,1:params.order+1)) ./ ...
                        (conj(bnAll(:,1:params.order+1)) .* bnAll(:,1:params.order+1) + params.regulConst));

    case 'softLimit'
        % Bernschuetz et al., Soft-Limiting der modalen Amplitudenverstärkung 
        gainLimLin = 10^(params.noiseGainDb / 20);
        radFilts = 2*gainLimLin/pi * abs(bnAll(:,1:params.order+1))./(bnAll(:,1:params.order+1)) ...
                    .* atan(pi ./ (2*gainLimLin*abs(bnAll(:,1:params.order+1))));

    case 'full'
        radFilts = 1 ./ bnAll(:,1:params.order+1);
        
    otherwise
        error('unkown radialFilter parameter');
        
end
