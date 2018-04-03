function [ thisval ] = obscenmoment(i, x)
%Compute i'th central moment of data array x.

xminusxmn = x - mean(x);
thisval = mean(xminusxmn.^i);

end

