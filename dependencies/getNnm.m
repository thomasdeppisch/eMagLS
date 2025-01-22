function Nnm = getNnm(N,zenRad,harmonicsDef)
% get the non-azimuthally dependent part of the SHs (SHs without CHs)
%
% Note that the Condon-Shortley phase of (-1)^m is not introduced in the
% code for the complex SH since it is included in the definition of the
% associated Legendre functions in Matlab (and it is canceled out in the code of the real SH).
arguments
    N (1,1)
    zenRad (1,1)
    harmonicsDef
end

Nnm = zeros((N+1)^2,1);
for nn=0:N
    Pn = legendre(nn,cos(zenRad)); % only gives results for positive m

    for mm=-nn:nn
        if strcmp(harmonicsDef,'complex')
            if mm<0
                Pnm = (-1)^abs(mm) * factorial(nn-abs(mm)) / factorial(nn+abs(mm)) * Pn(abs(mm)+1);
            else
                Pnm = Pn(mm+1);
            end
            
            Nnm(nn^2+nn+mm+1) = sqrt(((2*nn+1)*factorial(nn-mm))/((4*pi)*(factorial(nn+mm)))) * Pnm;

        elseif strcmp(harmonicsDef,'real')
            Nnm(nn^2+nn+mm+1) = (-1)^mm * sqrt(((2*nn+1)*factorial(nn-abs(mm)))/((4*pi)*(factorial(nn+abs(mm))))) * Pn(abs(mm)+1);
        end
    end
end