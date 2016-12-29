function y = dg(alpha,rm, r)
% y(t) = -alpha*(1+rm(t))^(-alpha-1)*r(t);
% size rm: T*1
% size r: T*k
T = length(rm);
y = -alpha.*ones(1,T)*(ones(T,1)+rm).^(-alpha-1)*r;