function [result] = cal_I(i,N)

if N == 1
    result = 0.5;
else
    if mod(i,2) == 1
        result = cal_I((i+1)/2,N/2)^2;
    else
        result = 2*cal_I(i/2,N/2) - cal_I(i/2,N/2)^2;
    end
end





