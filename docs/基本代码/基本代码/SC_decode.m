function [UR] = SC_decode(X,N,R,SNR)

%参数信息
%C : 接收码
%U : 信息集
%N : 码长
%R : 码率
%UL : 信息位长度
%SNR : 信噪比dB
%snr : 信噪比

snr = 10^(SNR/10);
UL = floor(N * R);
[num,~] = size(X);

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
UI = I(:,1:UL);
temp = zeros(1,N);
for i = 1 : UL
    temp(I(i)) = 1;
end

UR = zeros(num,UL);
for n = 1 : num
    %配置U
    U = zeros(1,N);

    %解码
    %判决第一个比特
    result = cal_L1(X(n,:),N,variance);
    u_e = judge(result);
    if temp(1) == 0
        U(1) = 0;
    else
        U(1) = u_e;
    end
        
    %判决剩余比特
    for i = 2 : N
        result = cal_LLR(X(n,:),U(1:i-1),N,i,variance);
        u_e = judge(result);
        if temp(i) == 0
            U(i) = 0;
        else
            U(i) = u_e;
        end
    end

    %恢复信息集
    for i = 1 : UL
        UR(n,i) = U(UI(i));
    end
end





