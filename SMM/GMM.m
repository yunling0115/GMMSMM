%% 1 GMM
clear all;
close;

% Moment equation: 0 = E[(c(t+1)/c(t))^(-gamma)*Re(t+1)]
% g = 1-E(mR)

%% Import data
    A = importfile('psgmmdata.xlsx', 'Sheet1', 'A2:F64');
    year = A(:,1);
    growth = A(:,2);
    mktrf = A(:,3);
    smb = A(:,4);
    hml = A(:,5);
    rf = A(:,6);

%% 1. Estimation
    
    n=5; % window size
    
    % (1) use 'rm-rf'
    exret1 = mktrf;
    % a. use identity matrix - first stage
    W1 = eye(1,1); 
    [b11,fval11] = fminsearch(@(gamma) gWg(growth,exret1,W1,gamma), 0); 
    S11 = S(growth, exret1, n, b11);
    % b. use S - second stage
    S11;
    [b12,fval12] = fminsearch(@(gamma) gWg(growth,exret1,S11,gamma), b11); 
    S12 = S(growth, exret1, n, b12);

    % (2) use 'rm-rf' and 'hml'
    exret2 = [mktrf, hml];
    % a. use identity matrix - first stage
    W2 = eye(2,2);
    [b21,fval21] = fminsearch(@(gamma) gWg(growth,exret2,W2,gamma), 0); 
    S21 = S(growth, exret2, n, b21);
    % b. use S - second stage
    S21;
    [b22,fval22] = fminsearch(@(gamma) gWg(growth,exret2,S21,gamma), b21); 
    S22 = S(growth, exret2, n, b22);

    % (3) use 'rm-rf' and 'smb'
    exret3 = [mktrf, smb];
    % a. use identity matrix - first stage
    W3 = eye(2,2);
    [b31,fval31] = fminsearch(@(gamma) gWg(growth,exret3,W3,gamma), 0); 
    S31 = S(growth, exret3, n, b31);
    % b. use S - second stage
    S31;
    [b32,fval32] = fminsearch(@(gamma) gWg(growth,exret3,S31,gamma), b31); 
    S32 = S(growth, exret3, n, b32);

%% 2. Standard errors
    
    delta = 1e-2; % finite difference
    
    std11 = SD(growth, exret1, W1, S11, b11, delta);
    std12 = SD(growth, exret1, S11, S12, b12, delta);
    
    std21 = SD(growth, exret2, W2, S21, b21, delta);
    std22 = SD(growth, exret2, S21, S22, b22, delta);
    
    std31 = SD(growth, exret3, W3, S31, b31, delta);
    std32 = SD(growth, exret3, S31, S32, b32, delta);   

%% Report gamma_hat, gamma_std_hat, and gWg
    
    clc;
    
    b11
    b12
    b21
    b22
    b31
    b32
    
    fval11
    fval12
    fval21
    fval22
    fval31
    fval32
    
    std11
    std12
    std21
    std22
    std31
    std32      

%% 3. Plot of g'Wg vs. gamma

    gamma= 0:0.1:100;
    for i=1:length(gamma);
        ga = gamma(i);
        gWg11(i) = gWg(growth,exret1,W1,ga);
        gWg12(i) = gWg(growth,exret1,S11,ga);
        gWg21(i) = gWg(growth,exret2,W2,ga);
        gWg22(i) = gWg(growth,exret2,S21,ga);
        gWg31(i) = gWg(growth,exret3,W3,ga);
        gWg32(i) = gWg(growth,exret3,S31,ga);
    end;
    
    subplot(1,2,1);  
    plot(gamma,gWg11,'b',gamma,gWg21,'r',gamma,gWg31,'g');
    title('Identity Matrix');
    xlabel('gamma');
    ylabel('gWg');
    subplot(1,2,2);  
    plot(gamma,gWg12,'b',gamma,gWg22,'r',gamma,gWg32,'g');
    title('Optimal Weighting Matrix');
    xlabel('gamma');
    ylabel('gWg');
    % mktrf
    figure(1);
    subplot(1,2,1); plot(gamma,gWg11,'b'); title('Identity Matrix');
    subplot(1,2,2); plot(gamma,gWg12,'b'); title('Optimal Weighting Matrix');   
    % mktrf + hml
    figure(2);
    subplot(1,2,1); plot(gamma,gWg21,'r'); title('Identity Matrix');
    subplot(1,2,2); plot(gamma,gWg22,'r'); title('Optimal Weighting Matrix');   
    % mktrf + smb
    figure(3);
    subplot(1,2,1); plot(gamma,gWg31,'g'); title('Identity Matrix');
    subplot(1,2,2); plot(gamma,gWg32,'g'); title('Optimal Weighting Matrix');   


%% 4. Plot of g vs. gamma

    gamma= 0:0.1:100;
    for i=1:length(gamma);
        ga = gamma(i);
        g11(i,:) = g(growth,exret1,W1,ga);
        g12(i,:) = g(growth,exret1,S11,ga);
        g21(i,:) = g(growth,exret2,W2,ga);
        g22(i,:) = g(growth,exret2,S21,ga);
        g31(i,:) = g(growth,exret3,W3,ga);
        g32(i,:) = g(growth,exret3,S31,ga);
    end;
    % mktrf
    figure(4);
    subplot(1,2,1); plot(gamma,g11,'b'); title('Identity Matrix');
    subplot(1,2,2); plot(gamma,g12,'b'); title('Optimal Weighting Matrix');   
    % mktrf + hml
    figure(5);
    subplot(1,2,1); plot(gamma,g21(:,1),'r',gamma,g21(:,2),'r-.'); title('Identity Matrix');
    subplot(1,2,2); plot(gamma,g22(:,1),'r',gamma,g22(:,2),'r-.'); title('Optimal Weighting Matrix');
    % mktrf + smb
    figure(6);
    subplot(1,2,1); plot(gamma,g31(:,1),'g',gamma,g31(:,2),'g-.'); title('Identity Matrix');
    subplot(1,2,2); plot(gamma,g32(:,1),'g',gamma,g32(:,2),'g-.'); title('Optimal Weighting Matrix');   

%% 5. Plot of min_g'Wg vs. min_gamma

    figure(7);
    plot(b11,fval11,'*'); hold on;
    plot(b12,fval12,'o');
    plot(b21,fval21,'r*');
    plot(b22,fval22,'ro');
    plot(b31,fval31,'g*');
    plot(b32,fval32,'go'); hold off;




