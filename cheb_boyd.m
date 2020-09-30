function [D1, D2, x] = cheb_boyd(N, pf)
% CHEB_BOYD - Compute differential matrix 
% for Chebyshev collocation method;
% It contains an optional parity factor (pf)
t = (pi/(2*N)*(1:2:(2*N-1)))';
x = cos(t); n = 0:(N-1);
ss = sin(t); cc = cos(t); 
sx = repmat(ss,1,N); cx = repmat(cc,1,N);
nx = repmat(n,N,1);  tx = repmat(t,1,N);
tn = cos(nx.*tx);
if  pf==0   
    phi2 = tn;              
    PT = -nx.*sin(nx.*tx);
    phiD2 = -PT./sx;        
    PTT = -nx.^2.*tn;
    phiDD2 = (sx.*PTT-cx.*PT)./sx.^3; 
else
    phi2 = tn.*sx;
    PT = -nx.*sin(nx.*tx).*sx + tn.*cx;
    phiD2 = -PT./sx;
    PTT = -nx.^2.*tn.*sx ...
          - 2*nx.*sin(nx.*tx).*cx - tn.*sx;
    phiDD2 = (sx.*PTT-cx.*PT)./sx.^3;
end
D1 = phiD2 /phi2; % the first derivatives
D2 = phiDD2/phi2; % the second derivatives