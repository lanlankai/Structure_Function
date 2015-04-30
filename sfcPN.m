% [sfP,sfN,sfO,sfM,Nk]=sfcPN(x,tlag,norder,dq);
% This function is to estimate the structure function with the term of
% positive and negative increments
% Input
% x is one dimensional time series to be analyzed
% tlag is time separations
% norder is the maximum statistical order (the length must be one)
% dq is the increment of statistical order
% Output
%  S_1(\tau) (sfP) is the structure function with positive contribution
%  S_2(\tau) (sfN) is the structure function with negative contribution
%  S_3(\tau) (sfO) is the structure function with all contribution (|P|+|N|)
%  S_4(\tau) (sfM) is the structure function with all contribution (|P|-|N|)
% Nk is the number of each part
%    Nk(1,:) is the number of positive part of each separation
%    Nk(2,:) is the number of negative part of each separation
%    Nk(3,:) is the total number of each separation
% 
% To see the result: loglog(tlag, sf(i,:))
% 
% Written by Yongxiang HUANG 28/03/2010
% 
% See also: sf