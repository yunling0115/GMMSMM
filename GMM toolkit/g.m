function y = g(alpha,rm, r)
% y(t) = (1+rm(t))^(-alpha)*r(t);
% size rm: T*1
% size r: T*k
T = length(rm);
y = ones(1,T)*(ones(T,1)+rm).^(-alpha)*r;
