function J = getChToShExpansionMatrix(order, harmonicsDef)
% build a matrix that expands a set of circular harmonics to equatorial
% harmonics
%
% see eq. (19) in Ahrens et al., "﻿Spherical harmonic decomposition of a 
% sound field based on observations along the equator of a rigid spherical 
% scatterer", JASA, 2021.
%
% td 2024

J = zeros((order+1)^2, 2*order+1);
Nnm = getNnm(order, pi/2, harmonicsDef);
for n = 0 : order
    for m = -n : n 
        acnIdx = n^2+n+m+1;
        J(n^2+n+m+1, 2*abs(m)-(sign(m)==-1)+1) = Nnm(acnIdx);
    end
end


end


