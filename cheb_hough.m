% CHEB_HOUGH - Compute Hough functions 
% using Chebyshev collocation methods
clear; format long e
a = 6.370d6; g = 9.81d0;
omega = 2.d0*pi/(24.d0*3600.d0);
%s = 1.d0; sigma = 0.4986348375d0; % DW1
 s = 1.d0; sigma = 0.5d0;    % DW1
%s = 2.d0; sigma = 1.0d0;    % SW2
%s = 3.d0; sigma = 1.5d0;    % TW3
parity_factor = mod(s,2);
N = 62; [D1,D2,x] = cheb_boyd(N,parity_factor);
a2 = (1-x.^2)./(sigma^2-x.^2);
a1 = 2.*x.*(1-sigma^2)./(sigma^2-x.^2).^2;
a0 = -1./(sigma^2-x.^2).*((s/sigma) ...
     .*(sigma^2+x.^2)./(sigma^2-x.^2) ...
     +s^2./(1-x.^2));
A = diag(a2)*D2 + diag(a1)*D1 + diag(a0); 
[v,d] = eig(A); lamb = real(diag(d));
% sort eigenvalues and -vectors
[foo,ii] = sort(-lamb); 
lamb = lamb(ii); hough = real(v(:,ii));
% equivalent depth (km)
h = -4.d0*a^2*omega^2/g./lamb/1000.d0;
% compute Hough functions for wind components
b1 = (sigma^2-x.^2).*sqrt(1.d0-x.^2);
b2 = sqrt(1.d0-x.^2)./(sigma^2-x.^2);
hough_u = diag(s./b1)*hough ...
          - diag(b2.*x./sigma)*D1*hough;         
hough_v = diag((s/sigma).*x./b1)*hough ...
          - diag(b2)*D1*hough;
clf % plot Hough functions
for j = 1:60
u = hough(:,j); subplot(10,6,j)
plot(x, u,'LineWidth',2), grid on
end