% NALP_HOUGH - Compute Hough functions 
% using normalized associated Legendre 
% polynomials (ALP)
clear; format long e
a = 6.370d6; g = 9.81d0;
omega = 2.d0*pi/(24.d0*3600.d0);
%s = 1.d0; sigma = 0.4986348375d0; % DW1
 s = 1.d0; sigma = 0.5d0;    % DW1
%s = 2.d0; sigma = 1.0d0;    % SW2
%s = 3.d0; sigma = 1.5d0;    % TW3
N = 62; N2 = N/2; sf = s/sigma;
% define L(r) and M(r)
L = zeros(N,1); M = zeros(N,1);
for r = s:N+s-1
i = r-s+1;
% define L(r)
L(i) = sqrt((r+s+1)*(r+s+2)*(r-s+1)*(r-s+2))...
       /((2*r+3)*sqrt((2*r+1)*(2*r+5))...
       *(sf-(r+1)*(r+2)));
% define M(r)
if (s == 2) && (r == 2)
   M(i) = -(sigma^2*(sf-r*(r+1)))...
          /((r*(r+1))^2)...
          +(r+2)^2*(r+s+1)*(r-s+1)...
          /((r+1)^2*(2*r+3)*(2*r+1)...
          *(sf-(r+1)*(r+2))); 
else
   M(i) = -(sigma^2*(sf-r*(r+1)))...
          /((r*(r+1))^2)...
          +(r+2)^2*(r+s+1)*(r-s+1)...
          /((r+1)^2*(2*r+3)*(2*r+1)...
          *(sf-(r+1)*(r+2)))...
          +(r-1)^2*(r^2-s^2)...
          /(r^2*(4*r^2-1)*(sf-r*(r-1)));
end % if
if (M(i) == inf), M(i) = realmax; end
end % for
% build F1 & F2 matix
f1 = zeros(N2,N2); f2 = zeros(N2,N2);
for i = 1:N2
f1(i,i) = M(2*i-1);
f2(i,i) = M(2*i);
if (i+1 <= N2)
   f1(i,i+1) = L(2*i-1);
   f1(i+1,i) = L(2*i-1);
   f2(i,i+1) = L(2*i);
   f2(i+1,i) = L(2*i);   
end % if
end % for
% symmetric modes
[v1,d1] = eig(f1); lamb1 = diag(d1);
[~,ii] = sort(-lamb1);                 
lamb1 = lamb1(ii); v1 = v1(:,ii);    
ht1 = 4.d0*a^2*omega^2/g.*lamb1/1000.d0;
% anti-symmetric modes
[v2,d2] = eig(f2); lamb2 = diag(d2);
[~,ii] = sort(-lamb2);                 
lamb2 = lamb2(ii); v2 = v2(:,ii);    
ht2 = 4.d0*a^2*omega^2/g.*lamb2/1000.d0;
% Legendre-Gauss quadrature points
nlat = 94; [x,w] = lgwt(nlat,-1,1);
% normalized associated Legendre functions
prs = pmn_polynomial_value(nlat,N+s,s,x);
% compute Hough modes
h1 = zeros(nlat,N2); h2 = zeros(nlat,N2);
for i = 1:N2
for j = 1:N2
i1 = 2*j+s-1; i2 = 2*j+s;
for ii = 1:nlat
% symmetric modes
h1(ii,i) = h1(ii,i) + v1(j,i)*prs(ii,i1);
% anti-symmetric modes
h2(ii,i) = h2(ii,i) + v2(j,i)*prs(ii,i2);
end
end
end
% put them together
lamb = zeros(N,1); hough = zeros(nlat,N);
for i = 1:N2
for j = 1:nlat    
i1 = 2*i-1; i2 = 2*i;
lamb(i1) = lamb1(i); 
lamb(i2) = lamb2(i);
hough(j,i1) = h1(j,i); 
hough(j,i2) = h2(j,i);
end
end
[~,ii] = sort(1./lamb);  
lamb = lamb(ii); hough = hough(:,ii);
% equivalent depth (km)
h = 4.d0*a^2*omega^2/g.*lamb/1000.d0;
% compute Hough functions for wind components
b1 = (sigma^2-x.^2).*sqrt(1.d0-x.^2);
b2 = sqrt(1.d0-x.^2)./(sigma^2-x.^2);
dhdx = central_diff(hough,x);
hough_u = diag(s./b1)*hough ...
          - diag(b2.*x./sigma)*dhdx;        
hough_v = diag((s/sigma).*x./b1)*hough ...
          - diag(b2)*dhdx;
clf % plot Hough functions
for j = 1:60
u = hough(:,j); subplot(10,6,j)
plot(x, u,'LineWidth',2), grid on
end  