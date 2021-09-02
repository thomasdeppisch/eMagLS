function Hdiff=sph_besselh_diff(nmax,kr)
% derivative of the spherical hankel function of the second kind

H=sph_besselh(nmax+1,kr);
Hdiff=zeros(size(H,1),nmax+1);
ofs=1;
for n=0:nmax
   Hdiff(:,n+ofs)=n*H(:,n+ofs)./kr(:) - H(:,n+1+ofs);
end