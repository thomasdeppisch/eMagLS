function Y = getCH(N, aziRad, basisType)
% This function returns complex or real-valued circular harmonics, up to
% maximum order N, evaluated at azimuthal directions in radians aziRad.
% Y .. CHs of size numDirs x 2N+1, CHs C_m are ordered as in ACN:
% [C_0, C_-1, C_1, C_-2, C_2, ..., C_-N, C_N].
%
arguments
    N (1,1)
    aziRad (:,1)
    basisType {mustBeMember(basisType,{'real','complex'})}
end

Ndirs = size(aziRad, 1);
Nharm = 2*N+1;
Y = zeros(Ndirs, Nharm);

Y(:,1) = 1;
for nn = 1:N
    chIdxNegM = 2*nn;
    chIdxPosM = 2*nn+1; 
    
    if isequal(basisType, 'real')
        Y(:,chIdxNegM) = sqrt(2) * sin(nn*aziRad);
        Y(:,chIdxPosM) = sqrt(2) * cos(nn*aziRad);
    elseif isequal(basisType, 'complex')
       Y(:,chIdxNegM) = exp(-1i*nn*aziRad);
       Y(:,chIdxPosM) = exp(1i*nn*aziRad);
    end
end