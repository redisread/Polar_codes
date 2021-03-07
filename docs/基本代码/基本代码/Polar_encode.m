function [X] = Polar_encode(U,N,SNR)

%参数信息
%U : 信息集
%N : 码长
%R : 码率
%UL : 信息位长度
%SNR : 信噪比dB
%snr : 信噪比

snr = 10^(SNR/10);
[num,UL] = size(U);

%高斯近似法挑选信息集
variance = 1/snr;
n = log2(N);
LLR = zeros(n+1,N);
LLR(1,:) = 2/variance;
for i = 2 : n+1
    for j = 1 : N
        if mod(j,2) == 1
            LLR(i,j) = SGA( LLR( i-1,(j+1)/2 ) );
        else
            LLR(i,j) = 2*LLR(i-1,j/2);
        end
    end
end

[~,I] = sort(LLR(n+1,:),2,'descend');
% x = 1 : 1 : N;
% plot(x,I,'.');
UI = I(:,1:UL);

%构建码字
X = zeros(num,N);
for n = 1 : num
    for i = 1 : UL
        X(n,UI(i)) = U(n,i);
    end
end

%计算生成矩阵GN
GN = cal_GN(N);

%编码
X = mod(X*GN,2);

%符号化及叠加噪声
X = awgn(1-2*X,SNR);
% C = 1-2*C;








