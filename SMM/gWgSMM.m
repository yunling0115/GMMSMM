function OBJ = gWgSMM(W,theta)
% calculate g for SMM
global delta lnz z Pi r Z T2 indexZ;
global step k_min k_max k;
global Moments;
qq = q(theta);
% g = [mean(qq)-Moments(1), std(qq)-Moments(2)];
g = [mean(qq)-Moments(1), mean(qq.^2)-sum(Moments.^2)];
OBJ = g*(W\g');