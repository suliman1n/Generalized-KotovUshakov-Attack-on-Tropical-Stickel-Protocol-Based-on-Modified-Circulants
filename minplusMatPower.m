function [C] = minplusMatPower(A,p)
%TropMatPower - Completes Tropical matrix multiplication p times (to the
%power of p)
%   Detailed explanation goes here
[n,m] = size(A);
if (n~=m)
    error("Dimension Error! Not a square matrix")
end 
if p ==0
    C = minplusTropId(n);
    return
end 
temp = A;
for i =2:p
    temp = minplusMulti(A,temp);
end
C = temp;
end