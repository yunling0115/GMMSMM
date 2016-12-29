function SS = S(growth, exret, n, b_hat)
% S = sum(j=-k:k) w(|j|)/T* sum(t=k+1:T-k)
% {[u(t,b1)-g][u(t-j,b1)-g]'}: N by N
[T,N] = size(exret);
wt = [0:n,n-1:-1:0]'; % 2k+1 by 1
% Or wt = bartlett(n+1)
u = (growth.^(-b_hat)*ones(1,N)).*exret; % T by N
g = 1/T.*sum(u,1); % 1 by N
u_dm = u-ones(T,1)*g; % demeaned u: T by N
wt_expand = [zeros(T-n-1,1);bartlett(2*n+1);zeros(T-n-1,1)];
SSr =  1/T*sum((wt_expand*ones(1,N^2)).*xcov(u_dm));
SS = reshape(SSr,[],N);



