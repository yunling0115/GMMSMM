function std = SD(growth, exret, W, SS, b_hat, delta)
% Var(b_hat) = 1/T*inv(d'Wd)*(d'WSWd)*inv(d'Wd)
% d = g' at b_hat
[T,N] = size(exret);
u1 = (growth.^(-b_hat)*ones(1,N)).*exret; % T by N
g1 = 1/T.*sum(u1,1); % 1 by N
b_hat = b_hat+delta;
u2 = (growth.^(-b_hat)*ones(1,N)).*exret; % T by N
g2 = 1/T.*sum(u2,1); % 1 by N
d = (g2-g1)/delta; % 1 by N
% Or use
% d1 = -(log(growth).*exp(-b_hat*log(growth)))*ones(1,N).*exret;
% d = 1/T.*sum(d1,1);
std = sqrt(inv(d*W*d')*(d*W*SS*W*d')*inv(d*W*d')/T);


