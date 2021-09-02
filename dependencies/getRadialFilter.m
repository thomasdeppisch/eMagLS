function radFilts = getRadialFilter(params)
% return radial filters
% Thomas Deppisch, 2021

c = 343;
nfft = params.oversamplingFactor*params.irLen;
f = linspace(0,params.fs/2,nfft/2+1)';
k = 2*pi*f/c;
kr = k * params.smaRadius;

if strcmp(params.waveModel, 'pointSource')
    krSource = k * params.sourceDist; % in this case params.sourceDist needs to be set by the user
end

% set radial filtering for waveModel + arrayType combination
switch params.arrayType
    case 'rigid'
        bn = @(N_,kr_) 1i ./ ((kr_.').^2 .* sph_besselh_diff(N_, kr_).');

    case 'open'
        bn = @(N_,kr_) sph_besselj(N_, kr_).';
        
    case 'directional' % open array with first-order directional mics -> see Politis, Array Response Simulator
        % 0 .. omni, 0.5 .. cardioid, 1 .. fig-of-eight
        bn = @(N_,kr_) (params.dirCoeff*sph_besselj(N_, kr_).' - 1i*(1-params.dirCoeff)*sph_besselj_diff(N_, kr_).');

    otherwise
        error('unkown arrayType parameter');
end

n = (0:params.order);
switch params.waveModel
    case 'planeWave'
        bnAll = 4*pi*1i.^n .* bn(params.order, kr).';

    case 'pointSource'
        hnAll = sph_besselh(params.order, krSource);
        bnAll = 4*pi*(-1i) .* k .* hnAll .* bn(params.order, kr).';

    otherwise
        error('unkown waveModel parameter');
end

    
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
