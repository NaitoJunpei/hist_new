function [rate_func] = HMM(x)
% [rate_func] = HMM(x,N)
%
% Function `HMM' returns optimal number of bins in a histogram
% used for density estimation.
% An assumption made is merely that samples are drawn from the density
% independently each other.
%
% Original paper:
% Mochizuki and Shinomoto, Analog and digital codes in the brain
% https://arxiv.org/abs/1311.4035
%
% Example usage:
% rate_func = HMM(x);
%
% Input argument
% x:    Sample data vector.
%
% Output argument
% rate_func: ratefunction.
%            2D array stores
%            1: begining time of each bins in second
%            2: rate of each bin
%
% Copyright (c) 2017, Kazuki Nakamura All rights reserved.
% shinomoto@scphys.kyoto-u.ac.jp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Setting
x = reshape(x,1,numel(x));
x_min = min(x);
x_max = max(x);

bin_width = ((x_max-x_min)/(length(x))) * 5;	% bin width = ISI * 5
EMloop_num=5000;        % number of EM itteration
mat_A=[0.999 0.001; 0.001 0.999];
vec_pi=[0.5,0.5];
mean_rate=length(x)/(x(length(x))-x(1));
vec_lambda=[(mean_rate*0.75)*bin_width (mean_rate*1.25)*bin_width];
vec_spkt=x-x(1);
vec_Xi = zeros(length(vec_spkt));
bin_num=ceil(vec_spkt(length(vec_spkt))/bin_width);
for i=1:length(vec_spkt)
    bin_id=fix(vec_spkt(i)/bin_width)+1;
    if bin_id<bin_num
        vec_Xi(bin_id)=vec_Xi(bin_id)+1;
    end
end

% HMM_E_step

mat_emission=zeros(length(vec_Xi),length(vec_lambda));
for n=1:length(vec_Xi)
    for i=1:length(vec_lambda)
        mat_emission(n,i)=vec_lambda(i)^vec_Xi(n)*exp(-1.0*vec_lambda(i))/factorial(vec_Xi(n));
    end
end


num_of_states=length(vec_pi);
num_of_obs=length(mat_emission);

alpha_0 = mat_emission(1:num_of_states,1) * vec_pi;
C_0 = sum(alpha_0);

vec_C = C_0;
alpha_0 = alpha_0 / C_0;
mat_alpha=alpha_0;

for n=2:num_of_obs
    alpha_n = [];
    for i=1:num_of_states
        sum_j=0.0;
        for j=1:num_of_states
            sum_j=sum_j+mat_alpha(j,n-1)*mat_A(j,i);
        end
        alpha_n_i=mat_emission(n,i)*sum_j;
        alpha_n=[alpha_n;alpha_n_i];
    end
    C_n = sum(alpha_n);
    vec_C=[vec_C,C_n];
    alpha_n = alpha_n/C_n;
    mat_alpha=[mat_alpha,alpha_n];
end

% n=N-1
for i=1:num_of_states
    mat_beta(num_of_obs-1,i)=1.0;
end

%n<N-1
for m=2:num_of_obs-1
    n=(num_of_obs-m);
    for i=1:num_of_states
        sum_j=0.0;
        for j=1:num_of_states
            sum_j= sum_j+mat_beta(n+1,j)*mat_emission(n+1,j)*mat_A(i,j);
        end
        mat_beta(n,i)=(sum_j/vec_C(n+1));
    end
end

% gamma matrix
mat_Gamma = mat_alpha.'*mat_beta.';
% Xi matrix
for m=2:num_of_obs-2
    for i=1:num_of_states
        for j=1:num_of_states
            mat_Xi(m,i,j)=(mat_alpha(i,m)*mat_emission(m+1,j)*mat_A(i,j)*mat_beta(m+1,j))/vec_C(m+1);
        end
    end
end

%%%%%%%%%%%%%%%

mat_A_old=mat_A;
vec_pi_old=vec_pi;
vec_lambda_old=vec_lambda;
loop=0;
flag=0;
while(loop<=EMloop_num || flag==0)
    
    num_of_states=length(mat_A);
    num_of_obs=length(vec_Xi);
    
    % maximize wrt pi vector
    pi_denom=sum(mat_Gamma(1,1:num_of_states));
    vec_pi_new=mat_Gamma(1,1:num_of_states)/pi_denom;
    
    % maximize wrt lambda vector
    vec_lambda_new=zeros(num_of_states);
    for k=1:num_of_states
        lambda_denom = sum(mat_Gamma(:,k));
        lambda_nume = sum(mat_Gamma(:,k))*vec_Xi(n);
        
        if lambda_nume==0.0
            vec_lambda_new(k)=0.0;
        else
            vec_lambda_new(k)=lambda_nume/lambda_denom;
        end
    end
    
    % maximize wrt A matirx
    mat_A_new=zeros(num_of_states, num_of_states);
    for j=1:num_of_states
        A_denome = sum(sum(mat_Xi(:,j,:)));
        for k=1:num_of_states
            A_nume=sum(mat_Xi(:,j,k));
            if(A_nume==0.0)
                mat_A_new(j,k)=0.0;
            else
                mat_A_new(j,k)=A_nume/A_denome;
            end
        end
    end
    
    vec_pi=vec_pi_new;
    vec_lambda=vec_lambda_new;
    mat_A=mat_A_new;
    
    num_state=length(vec_pi);
    sum_check = sum(abs(mat_A_old - mat_A)) + sum(abs(vec_pi_old-vec_pi)) + sum(abs(vec_lambda_old - vec_lambda));
    if sum_check/(1.0*num_state*(num_state+2))<1.0e-7
        flag = flag + 1;
    end
    mat_A_old=mat_A;
    vec_pi_old=vec_pi;
    vec_lambda_old=vec_lambda;
    
    % HMM_E_step
mat_emission=zeros(length(vec_Xi),length(vec_lambda));
for n=1:length(vec_Xi)
    for i=1:length(vec_lambda)
        mat_emission(n,i)=vec_lambda(i)^vec_Xi(n)*exp(-1.0*vec_lambda(i))/factorial(vec_Xi(n));
    end
end


num_of_states=length(vec_pi);
num_of_obs=length(mat_emission);

alpha_0 = mat_emission(1:num_of_states,1) * vec_pi;
C_0 = sum(alpha_0);

vec_C = C_0;
alpha_0 = alpha_0 / C_0;
mat_alpha=alpha_0;

for n=2:num_of_obs
    alpha_n = [];
    for i=1:num_of_states
        sum_j=0.0;
        for j=1:num_of_states
            sum_j=sum_j+mat_alpha(j,n-1)*mat_A(j,i);
        end
        alpha_n_i=mat_emission(n,i)*sum_j;
        alpha_n=[alpha_n;alpha_n_i];
    end
    C_n = sum(alpha_n);
    vec_C=[vec_C,C_n];
    alpha_n = alpha_n/C_n;
    mat_alpha=[mat_alpha,alpha_n];
end

% n=N-1
for i=1:num_of_states
    mat_beta(num_of_obs-1,i)=1.0;
end

%n<N-1
for m=2:num_of_obs-1
    n=(num_of_obs-m);
    for i=1:num_of_states
        sum_j=0.0;
        for j=1:num_of_states
            sum_j= sum_j+mat_beta(n+1,j)*mat_emission(n+1,j)*mat_A(i,j);
        end
        mat_beta(n,i)=(sum_j/vec_C(n+1));
    end
end

% gamma matrix
mat_Gamma = mat_alpha.'*mat_beta.';
% Xi matrix
for m=2:num_of_obs-2
    for i=1:num_of_states
        for j=1:num_of_states
            mat_Xi(m,i,j)=(mat_alpha(i,m)*mat_emission(m+1,j)*mat_A(i,j)*mat_beta(m+1,j))/vec_C(m+1);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%
    
    loop = loop+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HMM Viterbi Section
mat_emission=zeros(length(vec_Xi),length(vec_lambda));
for n=1:length(vec_Xi)
    for i=1:length(vec_lambda)
        mat_emission(n,i)=vec_lambda(i)^vec_Xi(n)*exp(-1.0*vec_lambda(i))/factorial(vec_Xi(n));
    end
end
num_of_states=length(mat_A);
num_of_obs=length(vec_Xi);

mat_hs_seq=zeros(num_of_states, num_of_obs);

%n=0
vec_logp_seq = zeros(num_of_states);
mat_hs_seq(:,1)=j;
for j=1:num_of_states
    mat_hs_seq(j,1)=j;
    vec_logp_seq(j) = log(vec_pi(j)*mat_emission(1,j))/log(10);
end
%n>0
for n=1:num_of_obs
    % copy the seq.
    mat_hs_seq_buf=mat_hs_seq;
    vec_logp_seq_buf=vec_logp_seq;
    
    %  nth node->j
    for j=1:num_of_states
        % n-1th node->j
        % compute logp for i->j trans
        for i=1:num_of_states
            vec_h_logprob_i(i)=vec_logp_seq(i)+log(mat_emission(n,j)*mat_A(i,j))/log(10);
        end
        % get max logp
        max_element=max(vec_h_logprob_i);
        for sea=1:length(vec_h_logprob_i)
            if vec_h_logprob_i(sea)==max_element
                max_pos=sea;
            end
        end
        vec_logp_seq_buf(j)=max_element;
        mat_hs_seq_buf(j,:)=mat_hs_seq(max_pos,:);
        mat_hs_seq_buf(j,n)=j;
        
    end
    % updata the seq.
    mat_hs_seq=mat_hs_seq_buf;
    vec_logp_seq=vec_logp_seq_buf;
end
max_element=max(vec_logp_seq);

for sea=1:length(vec_logp_seq)
    if vec_logp_seq(sea)==max_element
        max_pos=sea;
    end
end
vec_hidden=mat_hs_seq(max_pos,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rate_func = zeros(length(vec_Xi),2);

c_time=0.0;
for n=1:length(vec_Xi)
    state_id=vec_hidden(n);
    rate_func(n,1)=round(c_time*100)/100.0;
    rate_func(n,2)=round(vec_lambda(state_id)*100)/(bin_width*100.0);
    c_time = c_time + bin_width;
end
rate_func
