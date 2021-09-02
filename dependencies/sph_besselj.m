function H=sph_besselj(nmax,kr)
% spherical bessel function
% see e.g. Rafaely, Fundamentals of Spherical Array Processing, eq. 2.29

n=0:nmax;
[Nu,Kr]=meshgrid(n+1/2,kr);
H=sqrt(pi/2./Kr).*besselj(Nu,Kr);

