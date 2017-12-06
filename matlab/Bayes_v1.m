function [time, rate] = Bayes_v1(X)
%Program for estimating the instantaneous rate and regularity of occurrence
%
%Input
%X:             time series data
%Output
%time:          rate time
%rate:          firing rate
%
%Reference
%Estimating instantaneous irregularity of neuronal firing.
%T. Shimokawa and S. Shinomoto, Neural Computation (2009) 21:1931-1951. 
%
%version 0.1
%2010-3-31 Takeaki Shimokawa
%Copyright (c) 2010, Takeaki Shimokawa All rights reserved.
%%%%%%%%%%
% Bayes_v1.m by Kazuki Nakamura
%%%%%%%%%%

max_value=max(X);
min_value=min(X);

ISI=diff(X);
mu=length(X)/(max_value-min_value);
beta0=power(mu,-3);
beta=EMmethod(ISI,beta0);
[EL_N, VL_N, COVL_N]=KalmanFilter(ISI,beta);
rate=EL_N.';
time=(X(1:length(X)-1)+X(2:length(X)))/2;
end

%%%%
% EMmethod
% estimates a parameter beta with the EM algorithm.

% arguments:
% ISI: inter spike interval
% beta0: initial value of the parameter beta

% returns the estimated beta.
%%%%
function beta_new = EMmethod(ISI, beta0)
N=length(ISI);
beta=0;
beta_new=beta0;

for j=0:100
    beta=beta_new;
    [EL_N, VL_N, COVL_N]=KalmanFilter(ISI,beta);
    
    beta_new=0;
    t0=0;
    
    indexes = find(ISI(1:N-1));
    kalman1=EL_N;
    kalman2=VL_N;
    kalman3=COVL_N.';
    beta_new=sum(kalman2(indexes+1)+kalman2(indexes)-2*kalman3(indexes)+(kalman1(indexes+1)-kalman1(indexes)).^2./ISI(indexes).');
    t0=N-1-length(indexes);
    beta_new=(N-t0-1)/(2*beta_new);
end
end

%%%%
% KalmanFilter
% estimates the rate of event occurrence with the Kalman filtering.
% arguments:
% ISI: inter spike interval
% beta: the parameter
% returns the rate of event occurrence.
%%%%
function [EL_N, VL_N, COVL_N]=KalmanFilter(ISI,beta)
    N = length(ISI);
    IEL = N / sum(ISI);
    IVL = power(IEL / 3, 2);
    A = IEL - ISI(1) * IVL;
    EL = zeros(N,2);
    VL = zeros(N,2);

    EL_N = zeros(N,1);
    VL_N = zeros(N,1);
    COVL_N = zeros(N,1);

    EL0 = EL(1);
    EL1 = EL(2);
    VL0 = VL(1);
    VL1 = VL(2);

    EL0(1) = (A + sqrt(A * A + 4 * IVL)) / 2;
    EL0i1 = EL0(1);
    VL0(1) =  1 / (1 / IVL + 1 / power(EL0i1, 2));
    VL0i1 = VL0(1);
    
    % prediction and filtering
    for i =1:N-1
        EL1(i) = EL0i1;
        EL1i = EL0i1;
        VL1(i) = VL0i1 + ISI(i) / (2 * beta);
        VL1i = VL1(i);

        A = EL1i - ISI(i + 1) * VL1i;
        EL0(i + 1) =  (A + sqrt(A * A + 4 * VL1i)) / 2;
        EL0i1 = EL0(i+1);
        VL0(i + 1) =  1 / (1 / VL1i + 1 / power(EL0i1, 2));
        VL0i1 = VL0(i+1);
    end
    % smoothing
    % EL_N(N) = EL_Ni1 = EL0(N)
    % VL_N(N) = VL_Ni1 = VL0(N)
    EL_N(N) =  EL0i1;
    EL_Ni1 = EL_N(N);
    VL_N(N) =  VL0i1;
    VL_Ni1 = VL_N(N);
    
    for i = 1:N-1
        %i = N - 2 - i;
        H = VL0(i) / VL1(i);

        EL_N(i) = EL0(i) + H * (EL_Ni1 - EL1(i));
        EL_Ni1 = EL_N(i);
        VL_N(i) = VL0(i) + H * H * (VL_Ni1 - VL1(i));
        VL_Ni1 = VL_N(i);
        % COVL_N(i) = H * VL_N[i + 1]
    end
    COVL_N = VL0(1:N-1) ./ VL1(1:N-1) .* VL_N(2:N).';

end

