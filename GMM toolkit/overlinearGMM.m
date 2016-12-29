function varargout = overlinearGMM(A,B,eps)
% model: g(t) = A(t)+B(t)*phi
% Returns:
% varargout{1} = phi_est: dim(theta) by 1
% varargout{2} = W_est = std(g): dim(g) by dim(g)
% varargout{3} = g_est: dim(g) by 1
% varargout{4} = std(phi): dim(phi) by dim(phi)
% varargout{5} = Hansen's J statistic
% Inputs:
% size(A) = (T,dim(g),1)
% size(B) = (T,dim(g),dim(phi))
% Updating rule: phi=inv(B(T)'WB(T))*(-B(T)'WA(T)), g=A(T)+B(T)*phi, 
% W=inv(cov(g)), W0=I
mean_A = mean(A)'; % size(mean_A)=(dim(g),1)
temp_B = mean(B); mean_B(:,:) = temp_B(1,:,:); % size(mean_B)=(dim(g),dim(phi))
[T,dim_g,dim_phi]=size(B);
W = eye(dim_g); lagged_W = W+eye(dim_g);
phi = inv(mean_B'*W*mean_B)*(-mean_B'*W*mean_A);
lagged_phi = phi+ones(dim_phi,1);
while (norm(lagged_phi-phi)>=eps)
    laaged_W = W;
    lagged_phi = phi;
    for t=1:T
        tp_A(:,1) = A(t,:)';
        tp_B(:,:) = B(t,:,:);
        g(t,:) = tp_A+tp_B*phi; % size(g)=(T,dim(g),1)=(T,dim(g))
    end
    W = inv(cov(g)); % size(W)=(dim(g),dim(g))
    phi = inv(mean_B'*W*mean_B)*(-mean_B'*W*mean_A);
end
varargout{1} = phi; % size = (dim(theta),1)
varargout{2} = W; % size = (dim(g),dim(g))
varargout{3} = mean_A + mean_B*phi; % size = (dim(g),1)
g_est = varargout{3}; 
varargout{4} = inv(mean_B'*W*mean_B); % size = (dim(phi),dim(phi))
varargout{5} = T*(g_est'*W*g_est); % Hansen's J test statistic


