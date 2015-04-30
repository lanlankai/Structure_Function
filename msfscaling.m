function sf=msfscaling(x,y,qm,qn,tau)
% sf=msfscaling(x,y,qm,qn,tau)
% This fucntion is to estimate the mixed structure function by considering a sign
% power
% S_q(\tau)=<signp(\Delta x_{\tau}(t),qm)signp(\Delta y_{\tau}(t),qn)>
% signp(x,q)=|x|^2*sign(x)
% Input 
% x is the first time series will be analyzed
% y is the second time series will be analyzed
% qm is the statistical order for x
% qn is the statistical order for y
% tau is the maximum time delay
% Output
% sf.P=sfP is the positive contribution 
% sf.N=sfN is the negative contribution
% sf.O=sfO is all contribution |P|+|N|
% sf.M=sfM is contribution |P|-|N|
% sf.Nk=Nk is the ratio of positive and negative;
% sf.tau=tlag;
% 
% Written by Yongxiang Huang 29/03/2011
% 
% See aslo: sfcPNc
% 

if nargin~=5
    error('Five inputs are required!')
end

if length(x)~=length(y)
    error('The length of two time series should be the same!')
end

if length(tau)==1 && tau>10
    tau=[2:2:8 fix(10.^[1:0.1:log10(tau)])];
end

[sfP,sfN,sfO,sfM,Nk]=msfc(x,y,qm,qn,tau);

norder=[qm qn];

sf.P=[sfP];
sf.N=[sfN];
sf.O=[sfO];
sf.M=[sfM];
sf.Nk=Nk;
sf.q=norder;
sf.tau=tau;

    