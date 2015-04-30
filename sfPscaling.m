function sf=sfPscaling(x,tau,norder,dq)
% sf=sfscaling(x,tau,norder,dq)
% This fucntion is to estimate the structure function (two point average) by considering a sign
% power
% S_q(\tau)=<signp(\Delta u+{\tau}(t),q)>
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
% sf.Nk=Nk;
% sf.q=norder;
% sf.tau=tlag;
% 
% Written by Yongxiang Huang 02/02/2009
% 
% See aslo:sfc, sfcPN
% 

%   References:
%   HUANG Y., SCHMITT F.G., LU Z. LIU Y. Arbitrary order Hilbert spectral analysis 
%  for time series possessing scaling statistics: a comparison study
%  Physical Review E (submitted)

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


nDecade=ceil(log10(tau));% get the length of time delay in decade
tlag=[2 4 6 8 fix(10.^[1:0.1:nDecade])]; % the time delay
x=x-mean(x);
[sfP,sfN,sfO,sfM,Nk]=sfcPP(x,tlag,norder,dq);
norder=dq:dq:norder;
sf.P=sfP;
sf.N=sfN;
sf.O=sfO;
sf.M=sfM;
sf.Nk=Nk;
sf.q=norder;
sf.tau=tlag;

    