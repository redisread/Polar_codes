function [result] = cal_LLR(y,u,N,i,v)

%参数信息
%y : 接收码
%u : 判决集
%N、i : 标号
%v : 方差
if N == 1
    result = (2/v)*y(1);
else
    if i == 1
        result = cal_L1(y,N,v);
    else
        if mod(i,2) == 1
            y1 = y(1:N/2);
            y2 = y((N/2)+1:N);
            u1 = zeros(1,(i-1)/2);
            u2 = zeros(1,(i-1)/2);
            for j = 1 : (i-1)/2
                u1(j) = mod(u(2*j-1)+u(2*j),2);
                u2(j) = u(2*j);
            end
            result = L_f( cal_LLR(y1,u1,N/2,(i+1)/2,v) , cal_LLR(y2,u2,N/2,(i+1)/2,v) );
        else
            y1 = y(1:N/2);
            y2 = y((N/2)+1:N);
            u1 = zeros(1,(i-2)/2);
            u2 = zeros(1,(i-2)/2);
            for j = 1 : (i-2)/2
                u1(j) = mod(u(2*j-1)+u(2*j),2);
                u2(j) = u(2*j);
            end
            result = L_g( cal_LLR(y1,u1,N/2,i/2,v) , cal_LLR(y2,u2,N/2,i/2,v) , u(i-1) );
        end
    end
end











