function [wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, order, ...
    shDefinition, shFunction)
% [wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, order, ...
%     shDefinition, shFunction)
%
% calculates least-squares decoding filters
% see Schoerkhuber, Zaunschirm, Hoeldrich,
% "Binaural Rendering of Ambisonic Signals via Magnitude Least Squares"
%
% wLsL                   .. time-domain decoding filter for left ear
% wLsR                   .. time-domain decoding filter for right ear
% hL                     .. HRIR set for left ear (numSamples x numDirections)
% hR                     .. HRIR set for right ear (numSamples x numDirections)
% hrirGridAziRad         .. grid azimuth angles in radians of HRIR set (numDirections x 1)
% hrirGridZenRad         .. grid zenith angles in radians of HRIR set (numDirections x 1)
% order                  .. SH output order
% shDefinition           .. SH basis type according to utilized shFunction, default: 'real'
% shFunction             .. SH basis function (see testEMagLs.m for example), default: @getSH
%
% This software is licensed under a Non-Commercial Software License 
% (see https://github.com/thomasdeppisch/eMagLS/blob/master/LICENSE for full details).
%
% Thomas Deppisch, 2021

if nargin < 7; shFunction = @getSH; end
if nargin < 6 || isempty(shDefinition); shDefinition = 'real'; end

fprintf('with @%s("%s") ... ', func2str(shFunction), shDefinition);
Y = shFunction(order, [hrirGridAziRad, hrirGridZenRad], shDefinition);
pinvY = pinv(Y);

wLsL = hL * pinvY';
wLsR = hR * pinvY';

end
