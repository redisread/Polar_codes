function [result] = cal_L1(y,N,v)

%参数信息
%y : 接收码
%N : 标号
%v : 方差

if N == 1
    result = (2/v)*y(1);
else
    y1 = y(1,1:N/2);
    y2 = y(1,(N/2)+1:N);
    result = L_f( cal_L1(y1,N/2,v) , cal_L1(y2,N/2,v) );
%     result = ( cal_L1(y1,N/2,v) * cal_L1(y2,N/2,v) + 1 )/( cal_L1(y1,N/2,v) + cal_L1(y2,N/2,v) );
end