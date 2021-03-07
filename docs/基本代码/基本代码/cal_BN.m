function [BN] = cal_BN(N)

%BN = RN(I2?B(N/2))
if N == 2
    BN = eye(2);
else
    RN = zeros(N,N);
    for i = 1 : N
        if mod(i,2) == 1
            RN(i,(i+1)/2) = 1;
        else
            RN(i,(i+N)/2) = 1;
        end
    end
    BN = RN * kron(eye(2),cal_BN(N/2));
end




