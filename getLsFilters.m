function [wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, order)
% [wLsL, wLsR] = getLsFilters(hL, hR, hrirGridAziRad, hrirGridZenRad, order, fs)
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
% fs                     .. sampling frequency
%
% This software is licensed under a Non-Commercial Software License 
% (see https://github.com/thomasdeppisch/eMagLS/blob/master/LICENSE for full details).
%
% Thomas Deppisch, 2021

shDefinition = 'complex'; % real or complex

Y = getSH(order, [hrirGridAziRad,hrirGridZenRad], shDefinition);
pinvY = pinv(Y');

wLsL = hL * pinvY;
wLsR = hR * pinvY;

