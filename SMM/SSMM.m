function SS = SSMM(theta1,n)
% Calculate S given theta1
global delta lnz z Pi r Z T2 indexZ;
global step k_min k_max k;
global Moments;
qq = q(theta1);
u = [qq-Moments(1),qq.^2-sum(Moments.^2)];
[T,N] = size(u);
wt_expand = [zeros(T-n-1,1);bartlett(2*n+1);zeros(T-n-1,1)];
SSr = 1/T*sum((wt_expand*ones(1,N^2)).*xcov(u));
SS = reshape(SSr,[],N);