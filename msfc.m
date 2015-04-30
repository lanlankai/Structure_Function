% [sfP,sfN,sfO,sfM,Nk]=sfcPN(x,y,qm,qn,tau)
% This function is to estimate the mixed structure function between x and y
% with statistical order qm and qn and time delay tau
% Input
% x is the first time series to be analyzed
% y is the second time series to be analyzed
% qm is the statistical order for x
% qn is the statistical order for y
% tau is the corresponding time delay
% Output
% sfP is the positive contribution 
% sfN is the negative contribution
% sfO is all contribution |P|+|N|
% sfM is contribution |P|-|N|
% Nk  is the ratio of the positive and negative
% 
% see also:sfcaling sfcPN.c


