function sf=sfscaling(x,tau,norder,dq)
% sf=sfscaling(x,tau,norder,dq)
% This fucntion is to estimate the structure function by considering a sign
% power
% S_q(\tau)=<signp(\Delta u_{\tau}(t),q)>
% signp(x,q)=|x|^2*sign(x)
% Input 
% x is the time series will be analyzed
% tau is the maximum time delay
% norder is the maximum statistic order
% dq is the step of order
% Output
% sf.P=sfP is the positive contribution 
% sf.N=sfN is the negative contribution
% sf.O=sfO is all contribution |P|+|N|
% sf.M=sfM is contribution |P|-|N|
% sf.Nk=Nk; number of points of positive and negative value
% sf.Q=norder;
% sf.Tau=tlag;
% 
% Written by Yongxiang Huang 02/02/2009
% 
% See aslo:sfc, sfcPN
% 

%   References:
%   1) Huang Y., SCHMITT F.G., Lu Z., Fougairolles P., Gagne, Y. Liu Y., Second-Order Structure Function in fully developed turbulence
%   Physical Review E  82, 026319 􏰌2010􏰋
%   2) HUANG Y., SCHMITT F.G.,Hermand J.P. ,Gagne Y., LU Z. LIU Y. Arbitrary order Hilbert spectral analysis 
%  for time series possessing scaling statistics: a comparison study
%  Physical Review E 84, 016208 (2011)

if nargin<2
    error('You should give two parameters at least!')
end
if nargin==2
    norder=6;
    dq=0.5;
end
if nargin==3
    dq=0.5;
end
if size(x,1)>size(x,2)
    x=x';
end

if length(tau)==1
    nDecade=ceil(log10(tau));% get the length of time delay in decade
    tlag=[2 4 6 8 fix(10.^[1:0.1:nDecade])]; % the time delay
else
    tlag=tau;
end

if isa(x,'single')
    x=double(x);
end

[sfP,sfN,sfO,sfM,Nk]=sfcPN(x,tlag,norder,dq);

norder=0:dq:norder;
tmp=(Nk(1,:)-Nk(2,:))./(Nk(1,:)+Nk(2,:));
Nlag=length(tlag);

sf.P=[ones(1,Nlag);sfP];
sf.N=[ones(1,Nlag);sfN];
sf.O=[ones(1,Nlag);sfO];
sf.M=[tmp;sfM];
sf.Nk=Nk;
sf.Q=norder;
sf.Tau=tlag;

    
