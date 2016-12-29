function gg = g(growth, exret, W, gamma)
[T,N] = size(exret);
u = (growth.^(-gamma)*ones(1,N)).*exret; % T by N
gg = 1/T.*sum(u,1); % 1 by N
