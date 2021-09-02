function H=sph_besselh(nmax,kr)
% spherical hankel function of the second kind (due to the conj)
% (apply complex conjugation to get first kind)
% see e.g. Rafaely, Fundamentals of Spherical Array Processing, eq. 2.30

n=0:nmax;
[Nu,Kr]=meshgrid(n+1/2,kr);
H=sqrt(pi/2./Kr).*conj(besselh(Nu,Kr));

