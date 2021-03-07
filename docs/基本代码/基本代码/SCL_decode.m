function [UR] = SCL_decode(C,N,R,SNR,L)

%参数信息
%C : 接收码
%U : 信息集
%N : 码长
%R : 码率
%L : 保留路径数
%CL : 当前路径数
%UL : 信息位长度
%SNR : 信噪比dB
%snr : 信噪比

snr = 10^(SNR/10);
UL = floor(N * R);
[num,~] = size(C);

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
    %配置PM和U
    PM = zeros(1,2*L);
    for i = 1 : 2*L
        PM(i) = 100000;
    end
    P1 = PM;
    U = zeros(2*L,N);

    %解码
    %判决第一个比特
    U(1,1) = 0;
    U(1+L,1) = 1;
    result = cal_L1(C(n,:),N,variance);
    u_e = judge(result);
    if u_e == 0
        if temp(1) == 0
            U(1+L,1) = 0;
        end
        PM(1) = 0;
        PM(1+L) = abs(result);
    else
        if temp(1) == 0
            U(1+L,1) = 0;
            PM(1) = abs(result);
            PM(1+L) = abs(result);
        else
            PM(1) = abs(result);
            PM(1+L) = 0;
        end
    end
    [PM,I] = sort(PM,2,'ascend');
    CL = 1;
    if temp(1) == 1
        CL = 2;
        if CL > L
            CL = L;
        end
    end
    PM(CL+1:2*L) = P1(CL+1:2*L);
    a = U;
    for i = 1 : CL
        a(i,:) = U(I(i),:);
    end
    U = a;

    %判决剩余比特
    for i = 2 : N
        U(L+1:2*L,:) = U(1:L,:);
        for j = 1 : CL
            U(j,i) = 0;
            U(j+L,i) = 1;
            PM(j+L) = PM(j);

            result = cal_LLR(C(n,:),U(j,1:i-1),N,i,variance);
            u_e = judge(result);
            if u_e == 0
                if temp(i) == 0
                    U(j+L,i) = 0;
                end
                PM(j+L) = PM(j+L) + abs(result);
            else
                if temp(i) == 0
                    U(j+L,i) = 0;
                    PM(j) = PM(j) + abs(result);
                    PM(j+L) = PM(j+L) + abs(result);
                else
                    PM(j) = PM(j) + abs(result);
                end
            end
        end
        if temp(i) == 1
            CL = 2*CL;
            if CL > L
                CL = L;
            end
        else
            PM(CL+1:2*L) = P1(CL+1:2*L);
        end
        [PM,I] = sort(PM,2,'ascend');
        a = U;
        for s = 1 : CL
            a(s,:) = U(I(s),:);
        end
        U = a;
    end

    %恢复信息集
    for i = 1 : UL
        UR(n,i) = U(1,UI(i));
    end
end