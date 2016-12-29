function OBJ = gWg(growth, exret, W, gamma);
% growth = c(t+1)/c(t): T by 1
% exret: T by N
% W: N by N
% u(t+1) = (c(t+1)/c(t))^(-gamma)*Re(t+1)
[T,N] = size(exret);
u = (growth.^(-gamma)*ones(1,N)).*exret; % T by N
g = 1/T.*sum(u,1);
OBJ = g*W*g';