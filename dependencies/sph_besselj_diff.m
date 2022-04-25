function Jdiff = sph_besselj_diff(nmax,kr)
% derivative of the spherical bessel function

J=sph_besselj(nmax+1,kr);
Jdiff=zeros(size(J,1),nmax+1);
ofs=1;
for n=0:nmax
   Jdiff(:,n+ofs)=n*J(:,n+ofs)./kr(:) - J(:,n+1+ofs);
end

end
