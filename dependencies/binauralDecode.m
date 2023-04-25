function binauralOut = binauralDecode(in, inFs, decodingFilterLeft, decodingFilterRight, decodingFilterFs, ...
    compensateDelay, signal, signalFs, horRotAngleRad)
    % decode SH-domain (Ambisonic) signal to binaural, optionally convolve 
    % with mono signal (in case ambisonic signal is a rir)
    % 
    % Thomas Deppisch, 2021
    
    %% some resampling
    if (nargin > 6)
        if (signalFs ~= inFs)
            disp('binauralDecode: resampling signal');
            signal = resample(signal, inFs, signalFs);
        end
    end
    
    if (decodingFilterFs ~= inFs)
        disp('binauralDecode: resampling decoding filter');
        decodingFilterLeft = resample(decodingFilterLeft, inFs, decodingFilterFs);
        decodingFilterRight = resample(decodingFilterRight, inFs, decodingFilterFs);
    end

    %% rotate horizontally?
    if (nargin > 8 && ~isempty(horRotAngleRad) && horRotAngleRad ~= 0)
        % TODO: This function requires 
        %       https://github.com/polarch/Higher-Order-Ambisonics.git
        in = rotateHOA_N3D(in, rad2deg(horRotAngleRad), 0, 0);
    end
    
    %% ambi binaural decoding
    numSamplesIn = size(in, 1);
    numHarmonics = size(in, 2);
    
    leftEarSig = zeros(numSamplesIn, 1);
    rightEarSig = zeros(numSamplesIn, 1);
    
    for ii = 1:numHarmonics
        leftEarSig = leftEarSig + fftfilt(decodingFilterLeft(:, ii), in(:, ii));
        rightEarSig = rightEarSig + fftfilt(decodingFilterRight(:, ii), in(:, ii));
    end
    
    %% optionally convolve with signal
    if (nargin > 6 && ~isempty(signal))
        leftEarSig = fftfilt(leftEarSig, signal(:,1));
        rightEarSig = fftfilt(rightEarSig, signal(:,1));
    end
    
    %% output
    binauralOut = [leftEarSig, rightEarSig];
    
    if nargin > 5 && compensateDelay
        % assume half filter length as delay
        del = size(decodingFilterLeft,1) / 2;
        binauralOut = binauralOut(del:end,:);
    end

    % force real output (may be relevant in case of complex SHs)
    if ~isreal(binauralOut)
        warning('discarding imaginary part with sum of [%.2g, %.2g] in rendering result.', ...
            sum(abs(imag(binauralOut(:,1)))), sum(abs(imag(binauralOut(:,2)))));
        binauralOut = real(binauralOut);
    end
end
