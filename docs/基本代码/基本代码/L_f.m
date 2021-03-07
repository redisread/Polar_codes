function [y] = L_f(a,b)

% y = (a*b + 1)/(a + b);%LR
% y = log( (1+exp(a+b)) / (exp(a)+exp(b)) );%LLR

%硬件友好性LLR
y = sign(a)*sign(b)*min([abs(a),abs(b)]);

