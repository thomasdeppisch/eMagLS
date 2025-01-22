function [wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, order, ...
    shDefinition, shFunction)
% [wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, order, ...
%     shDefinition, shFunction)
%
% This function calculates least-squares decoding filters for head related impulse response data sets.
% For more information, please refer to
%   Schörkhuber, Zaunschirm, and Hoeldrich,
%   “Binaural Rendering of Ambisonic Signals via Magnitude Least Squares,”
%   in Fortschritte der Akustik -- DAGA 2018, 2018, pp. 339–342.
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
% (see https://github.com/thomasdeppisch/eMagLS/blob/main/LICENSE for full details).
%
% Thomas Deppisch, 2023

if nargin < 7; shFunction = @getSH; end
if nargin < 6 || isempty(shDefinition); shDefinition = 'real'; end

Y_conj = shFunction(order, [hrirGridAziRad, hrirGridZenRad], shDefinition)';
Y_pinv = pinv(Y_conj);

wLsL = hL * Y_pinv;
wLsR = hR * Y_pinv;

end
