function q = q(theta)
% solving the moments
% Mom: N by 1 (here N = 2)
% 1. Given a value of theta, solve the firm problem using discrete state space dynamic programming
% methods
    global delta lnz z Pi r Z T1 T2 indexZ;
    global step k_min k_max k;
    v1 = zeros(size(k));
    v2 = zeros(size(k));
    v = [v1;v2];
    r = 1/0.96-1;
    iter = 1;
    % felicity function
    R1 = (z(1)*k.^theta+(1-delta)*k)*ones(size(k'))-ones(size(k))*k';
    R2 = (z(2)*k.^theta+(1-delta)*k)*ones(size(k'))-ones(size(k))*k';
    R1(find(R1<0)) = -inf; 
    R2(find(R2<0)) = -inf;
    [v_next,index] = max([R1;R2] + 1/(1+r)*kron(Pi,ones(size(k)))*[v1';v2'],[],2); % along the second dim: select max column
    while norm(v_next-v)>eps & iter<10000
        iter = iter+1;
        v = v_next;
        vr = reshape(v,[],2);
        [v_next,index] = max([R1;R2] + 1/(1+r)*kron(Pi,ones(size(k)))*vr',[],2);
    end
    indexr = reshape(index,[],2);
    i1 = k(indexr(:,1));
    i2 = k(indexr(:,2));
 % 2. Given Simulated Z, get K
    K = zeros(T2,1);
    V = zeros(T2,1);
    indexK = zeros(T2,1);
    indexK(1) = (floor(length(k)/2));
    K(1) = k(indexK(1));
    vr = reshape(v,[],2);
    % Get simulated K
    for t=2:T2;
         indexK(t) = indexr(indexK(t-1),indexZ(t));
         V(t) = vr(indexK(t),indexZ(t));
    end;
    K = k(indexK);
% 3. Calculate q = V/k (T by 1)
    q_all = V./K;
    q = q_all(T1+1:T2);
    
    
    
