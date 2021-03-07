clc;
clear;

N = 1024;
R = 1/2;
SNR = 2;

U = randi([0,1],1,N*R);
X = Polar_encode(U,N,SNR);
tic;
UR = SC_decode(X,N,R,SNR);
toc
num = 0;
for i = 1 : N*R
    if U(i) ~= UR(i)
        num = num + 1;
    end
end