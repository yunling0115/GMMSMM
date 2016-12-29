%% GMM and SMM
clear all;
close;

%% 1 GMM
% Moment equation: 0 = E[(c(t+1)/c(t))^(-gamma)*Re(t+1)]
% g = 1-E(mR)

% Import data
    A = importfile('psgmmdata.xlsx', 'Sheet1', 'A2:F64');
    year = A(:,1);
    growth = A(:,2);
    mktrf = A(:,3);
    smb = A(:,4);
    hml = A(:,5);
    rf = A(:,6);

% 1. Estimation

    % (1) use 'rm-rf'
    exret1 = mktrf;
    % a. use identitfy matrix
    W1 = eye(1,1); 
    [b11,fval] = fminsearch(@(gamma) gWg(growth,exret1,W1,gamma), 0); % b11 = 52.5236
    % b. use optimal matrix - exactly identified (same)
    k = 5;
    S1 = S(growth, exret1, k, b11);
    [b12,fval] = fminsearch(@(gamma) gWg(growth,exret1,S1,gamma), b11); % b12 = 52.5236

    % (2) use 'rm-rf' and 'hml'
    exret2 = [mktrf, hml];
    % use identitfy matrix
    W2 = eye(2,2);
    [b21,fval] = fminsearch(@(gamma) gWg(growth,exret2,W2,gamma), 0); % b21 = 56.3383
    % b. use optimal matrix
    k = 5;
    S2 = S(growth, exret2, k, b21);
    [b22,fval] = fminsearch(@(gamma) gWg(growth,exret2,S2,gamma), b21); % b22 = 47.5716

    % (3) use 'rm-rf' and 'smb'
    exret3 = [mktrf, smb];
    % use identitfy matrix
    W3 = eye(2,2);
    [b31,fval] = fminsearch(@(gamma) gWg(growth,exret3,W3,gamma), 0); % b31 = 52.7749
    % b. use optimal matrix
    k = 5;
    S3 = S(growth, exret3, k, b31);
    [b32,fval] = fminsearch(@(gamma) gWg(growth,exret3,S3,gamma), b31); % b32 = 60.2030
    
    clc;
    
    b11
    b12
    b21
    b22
    b31
    b32
    % b11 = 52.5236
    % b12 = 52.5236
    % b21 = 56.3383
    % b22 = 61.6433
    % b31 = 52.7749
    % b32 = 60.2030
   
% 2. Plot of g'Wg vs. gamma

    gamma= 0:0.1:100;
    for i=1:length(gamma);
        ga = gamma(i);
        gWg11(i) = gWg(growth,exret1,W1,ga);
        gWg12(i) = gWg(growth,exret1,S1,ga);
        gWg21(i) = gWg(growth,exret2,W2,ga);
        gWg22(i) = gWg(growth,exret2,S2,ga);
        gWg31(i) = gWg(growth,exret3,W3,ga);
        gWg32(i) = gWg(growth,exret3,S3,ga);
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


% 3. Plot of g vs. gamma

    gamma= 0:0.1:100;
    for i=1:length(gamma);
        ga = gamma(i);
        g11(i,:) = g(growth,exret1,W1,ga);
        g12(i,:) = g(growth,exret1,S1,ga);
        g21(i,:) = g(growth,exret2,W2,ga);
        g22(i,:) = g(growth,exret2,S2,ga);
        g31(i,:) = g(growth,exret3,W3,ga);
        g32(i,:) = g(growth,exret3,S3,ga);
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

% 4. Standard errors
    
    delta = 1e-2;
    std11 = SD(growth, exret1, W1, b11, delta);
    std12 = SD(growth, exret1, S1, b12, delta);
    std21 = SD(growth, exret2, W2, b21, delta);
    std22 = SD(growth, exret2, S2, b22, delta);
    std31 = SD(growth, exret3, W3, b31, delta);
    std32 = SD(growth, exret3, S3, b32, delta);
    std11
    std12
    std21
    std22
    std31
    std32
    % std11 = 126.0887
    % std12 = 44.4477
    % std21 = 112.1873
    % std22 = 41.3511
    % std31 = 126.0200
    % std32 = 38.6894  

%% 2 SMM

% Import Data

% 1. Solve the firm problem using discrete state space dynamic programming
% methods

% 2. Report the average investment rate (investment/capital stock), and its
% standard deviation 

% 3. Report the average Tobin's q, V/k and its standard deviation

% 4. Calculate these same moments using Compustat data, use the definition
% of Tobin's Q from Eisfeldt and Rampini (2006 JME) (exclude firms with
% less than 100m)

% 5. Using a grid search, find the value of theta that sets these data
% moment sas close as possible to the actual moments. Use both identity
% matrix and the optiaml weighting matrix. Discuss the difference between
% the two estiamtes.



