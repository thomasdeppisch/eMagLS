function fade_win = getFadeWindow(irLen, relFadeLen)
% return combined fade-in and fade-out window
% Hannes Helmholz, 2022

    if nargin < 2; relFadeLen = 0.15; end % relative length of result fading window
    
    n_fadein = round(relFadeLen * irLen);
    n_fadeout = round(relFadeLen * irLen);
    hannin = hann(2*n_fadein);
    hannout = hann(2*n_fadeout);
    fade_win = [hannin(1:end/2); ...
        ones(irLen-(n_fadein+n_fadeout),1); hannout(end/2+1:end)];

end
