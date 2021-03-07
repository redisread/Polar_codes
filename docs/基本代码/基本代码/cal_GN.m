function [GN] = cal_GN(N)

%GN = BN * F^logN
%计算BN
BN = cal_BN(N);

%计算GN
n = log2(N);
f = [1,0;1,1];
F = f;
for i = 1 : n-1
    F = kron(F,f);
end
GN = BN * F;








