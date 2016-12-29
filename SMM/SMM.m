%% 2 SMM

%% 1. Solve the firm problem using discrete state space dynamic programming
% methods

clear all;
close;

    global delta lnz z Pi r;
    theta = 0.6;
    delta = 0.15;
    lnz = [-0.2,0.2];
    z = exp(lnz);
    Pi = [0.6 0.4;0.4 0.6];
    % discretization of capital
    global step k_min k_max k;
    step = 0.1;
    k_min = 0.1;
    k_max = 50;
    k = [k_min:step:k_max]';
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
    [v_next,indexZ] = max([R1;R2] + 1/(1+r)*kron(Pi,ones(size(k)))*[v1';v2'],[],2); % along the second dim: select max column
    while norm(v_next-v)>eps & iter<10000
        iter = iter+1;
        v = v_next;
        vr = reshape(v,[],2);
        [v_next,indexZ] = max([R1;R2] + 1/(1+r)*kron(Pi,ones(size(k)))*vr',[],2);
    end
    indexr = reshape(indexZ,[],2);
    i1 = k(indexr(:,1));
    i2 = k(indexr(:,2));
    % Stationary distribution
    Pi_k1 = zeros(length(k));
    Pi_k2 = zeros(length(k));
    Pi_k1(indexr(:,1)+length(k)*[0:length(k)-1]') = 1; Pi_k1 = Pi_k1';
    Pi_k2(indexr(:,2)+length(k)*[0:length(k)-1]') = 1; Pi_k2 = Pi_k2';
    Pi = [kron(Pi(1,:),Pi_k1);kron(Pi(2,:),Pi_k2)]; % big transition matrix
    Pii = sparse(Pi);
    [eg,ev] = eigs(Pii'); 
    Pz = eg(:,1)/sum(eg(:,1)); % stationary probability (eigs: only look at nonzero)
 
    
% 2. Report the average investment rate (investment/capital stock), and
% its standard deviation
     
    invrate1 = (i1-(1-delta)*k)./k;
    invrate2 = (i2-(1-delta)*k)./k;
    % average investment rate
    Einv = Pz'*[invrate1;invrate2]; % 
    Pzr = reshape(Pz,[],2);
    % average investment rate (by state)
    Einv1 = (Pzr(:,1))'*invrate1/sum(Pzr(:,1)); % 
    Einv2 = (Pzr(:,2))'*invrate2/sum(Pzr(:,2)); % 
    % standard deviation
    Stdinv = sqrt(Pz'*[invrate1.^2;invrate2.^2]-Einv^2);  
    
% 3. Report the average Tobin's q, V/k and its standard deviation
    q = v./[k;k];
    Eq = Pz'*q;
    Stdq = sqrt(Pz'*(q.^2)-Eq^2);
    
    clc;
    Einv
    Einv1
    Einv2
    Stdinv    
    Eq
    Stdq
    
    % Einv = 0.1578
    % Einv1 = 0.0784
    % Einv2 = 0.2373
    % Stdinv = 0.1260
    % Eq = 4.2613
    % Stdq = 0.3108
    
%% 4. Calculate these same moments using Compustat data, use the definition
% of Tobin's Q from Eisfeldt and Rampini (2006 JME) (exclude firms with
% less than 100m)
    
clear all; % redefine to drop off the sparse matrix
close;
       
    global delta lnz z Pi r;
    theta = 0.6;
    delta = 0.15;
    lnz = [-0.2,0.2];
    z = exp(lnz);
    Pi = [0.6 0.4;0.4 0.6];
    % discretization of capital
    global step k_min k_max k;
    step = 0.1;
    k_min = 0.1;
    k_max = 50;
    k = [k_min:step:k_max]';
    v1 = zeros(size(k));
    v2 = zeros(size(k));
    v = [v1;v2];
    r = 1/0.96-1;
    
    % real q = MV/BV
    % BV = AT
    % MV = AT + CEQL - UCEQ  - TXDB
    % screening: delete if AT<0, or UCEQ<0, or TXDB<0,
    % if TXDB=. then TXDB=0
    % STD(q) (wt = MV) = mean(q in [0,5], wt)
    % (did in stata)
    % mean = .9197415  std =  .0565921 
    
    global Moments; 
    Moments = [2.1,1.1];
    
%% 5. Using a grid search, find the value of theta that sets these data
% moments as close as possible to the actual moments. Use both identity
% matrix and the optiaml weighting matrix. Discuss the difference between
% the two estiamtes.

    global T1 T2 indexZ Z;
    T2 = 1250;
    T1 = 250;
    % 0.simulate fundemantal shock realization
    jump = (rand(T2,1)>0.6);
    indexZ = zeros(T2,1); % = 1 or 2 for z state
    indexZ(1) = 1;
    for t = 2:T2;
        if jump(t)==0 indexZ(t)=indexZ(t-1);
        else if indexZ(t-1)==1 indexZ(t)=2;
            else indexZ(t)=1;
            end;
        end;
    end;   
    Z = z(indexZ);    
    % 1.for a given theta for model, simulate data from model, and compute moments
    q0 = q(theta);    
    % 2.Search on a grid of theta to find min_theta
    global grid_theta;
    grid_theta = 0:0.01:1;
    
    % a. identity matrix
    W = eye(2,2);
    for i=1:length(grid_theta);
        gWg1(i) = gWgSMM(W,grid_theta(i));
    end;
    [val1,inx1] = min(gWg1);
    thetaopt1 = grid_theta(inx1); 
    figure(1);
    plot(grid_theta, gWg1); hold on;
    plot(grid_theta(inx1),gWg1(inx1),'r*'); hold off;
    
    % b. first-stage matrix of thetaopt1
    n = 60; % window size
    S1 = SSMM(thetaopt1,n); % use first-stage estimated theta
    for i=1:length(grid_theta);
        gWg2(i) = gWgSMM(S1,grid_theta(i));
    end;
    [val2,inx2] = min(gWg2);
    thetaopt2 = grid_theta(inx2);
    figure(2);
    plot(grid_theta, gWg2); hold on;
    plot(grid_theta(inx2),gWg2(inx2),'r*'); hold off;
    S2 = SSMM(thetaopt2,n); 
        
    % c. first-stage matrix of grid_theta
    n = 60; % window size
    for i=1:length(grid_theta);
        gWg21(i) = gWgSMM(W,grid_theta(i));
        S21 = SSMM(grid_theta(i),n);
        gWg22(i) = gWgSMM(S1,grid_theta(i));
    end;
    [val22,inx22] = min(gWg22);
    thetaopt22 = grid_theta(inx22);
    figure(2);
    plot(grid_theta, gWg22); hold on;
    plot(grid_theta(inx22),gWg22(inx22),'r*'); hold off;
    S22 = SSMM(thetaopt22,n); 
    
    clc
    thetaopt1
    thetaopt2
    thetaopt22
    S1
    S2
    S22
















