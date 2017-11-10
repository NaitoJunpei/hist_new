function [rate_func] = HMM_v2(x)
% [rate_func] = HMM(x)
%
% Function `HMM' returns the firing rate selected as an alternative hidden state.
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
% made by Yasuhiro Mochizuki
% revised by Kazuki Nakamura
% Contact: Shigeru Shinomoto: shinomoto@scphys.kyoto-u.ac.jp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

onset = x(1) - 0.001 * (x(length(x)) - x(1));
offset = x(length(x)) + 0.001 * (x(length(x)) - x(1));
optw = (offset-onset)/(length(x)) * 5;
rate_func = get_hmm_ratefunc(x, optw);

drawHMM(rate_func)

end
%sub function

function vec_Xi = get_vec_Xi(vec_spkt, bin_width)
%%%%vec_Xi = zeros(length(vec_spkt),1);
bin_num=ceil(vec_spkt(length(vec_spkt))/bin_width);
vec_Xi = zeros(bin_num, 1);
for i=1:length(vec_spkt)
    bin_id=fix(vec_spkt(i)/bin_width)+1;
    if bin_id<bin_num + 1
        vec_Xi(bin_id) = vec_Xi(bin_id)+1;
    end
end
end

% func to get emisson probs.
function mat_emission = get_mat_emission(vec_Xi, vec_lambda)
mat_emission=zeros(length(vec_Xi),length(vec_lambda));
for n=1:length(vec_Xi)
    for i=1:length(vec_lambda)
        mat_emission(n,i)=vec_lambda(i)^vec_Xi(n)*exp(-1.0*vec_lambda(i))/factorial(vec_Xi(n));
    end
end
end

% func to get alpha and C
function [vec_C, mat_alpha] = get_alpha_C(mat_A, vec_pi, mat_emission)
num_of_states=length(vec_pi);
num_of_obs=length(mat_emission(:, 1));

vec_C = zeros(num_of_obs,1);
%n=1
for i=1:num_of_states
    mat_alpha(1,i) = mat_emission(1,i) * vec_pi(i);
end
vec_C(1) = sum(mat_alpha(1,:));
for i=1:num_of_states
    mat_alpha(1,i) = mat_alpha(1,i)/vec_C(1);
end

%n>1
for n=2:num_of_obs
    for i=1:num_of_states
        sum_j=0.0;
        for j=1:num_of_states
            sum_j = sum_j + mat_alpha(n-1,j)*mat_A(j,i);
        end
        mat_alpha(n,i) = mat_emission(n,i)*sum_j;
    end
    vec_C(n)=sum(mat_alpha(n,:));
    mat_alpha(n,:) = mat_alpha(n,:)./vec_C(n);
end
end

function mat_beta = get_beta(mat_A, vec_pi, mat_emission, vec_C)
num_of_states=length(vec_pi);
num_of_obs=length(mat_emission(:, 1));

% initialize
mat_beta = zeros(num_of_obs, num_of_states);

% n=N
mat_beta(num_of_obs,:)=1.0;

% n<N
for m=2:num_of_obs
    n=(num_of_obs+1-m);
    for i=1:num_of_states
        sum_j=0.0;
        for j=1:num_of_states
            sum_j= sum_j+mat_beta(n+1,j)*mat_emission(n+1,j)*mat_A(i,j);
        end
        mat_beta(n,i)=(sum_j/vec_C(n+1));
    end
end
end


function [mat_Gamma, mat_Xi] = get_Gamma_Xi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C)
[num_of_obs, num_of_states] = size(mat_emission);
% gamma matrix
mat_Gamma=zeros(num_of_obs, num_of_states);
for n=1:num_of_obs
    for i=1:num_of_states
        mat_Gamma(n,i)=mat_alpha(n,i)*mat_beta(n,i);
    end
end

% Xi matrix
mat_Xi=zeros(num_of_obs-1, num_of_states, num_of_states);
for m=1:num_of_obs-1
    for i=1:num_of_states
        for j=1:num_of_states
            mat_Xi(m,i,j)=(mat_alpha(m,i)*mat_emission(m+1,j)*mat_A(i,j)*mat_beta(m+1,j))/vec_C(m+1);
        end
    end
end
end

% HMM_E_step
function [mat_Gamma, mat_Xi] = HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi)
mat_emission=get_mat_emission(vec_Xi,vec_lambda);

[vec_C, mat_alpha]=get_alpha_C(mat_A, vec_pi, mat_emission);

mat_beta=get_beta(mat_A, vec_pi, mat_emission, vec_C);

[mat_Gamma, mat_Xi]=get_Gamma_Xi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C);
end

% HMM_M_step
function [vec_pi_new, vec_lambda_new, mat_A_new] = HMM_M_step(vec_Xi,  mat_A,  vec_lambda,  vec_pi,  mat_Gamma,  mat_Xi)
num_of_states=length(mat_A);
num_of_obs=length(vec_Xi);

% maximize wrt pi vector
pi_denom=sum(mat_Gamma(1,:));
vec_pi_new=mat_Gamma(1,:)./pi_denom;

% maximize wrt lambda vector
vec_lambda_new=zeros(num_of_states,1);
for k=1:num_of_states
    lambda_denom=sum(mat_Gamma(:,k));
    
    lambda_nume=0.0;
    for n=1:num_of_obs
        lambda_nume=lambda_nume+(mat_Gamma(n,k)*vec_Xi(n));
    end
    
    if lambda_nume==0.0
        vec_lambda_new(k)=0.0;
    else
        vec_lambda_new(k)=lambda_nume/lambda_denom;
    end
end

% maximize wrt A matirx
mat_A_new=zeros(num_of_states, num_of_states);
for j=1:num_of_states
    A_denome=0.0;
    for n=1:num_of_obs-1
        for l=1:num_of_states
            A_denome=A_denome+mat_Xi(n,j,l);
        end
    end
    for k=1:num_of_states
        A_nume=0.0;
        for n=1:num_of_obs-1
            A_nume = A_nume+mat_Xi(n,j,k);
            if A_nume==0.0
                mat_A_new(j,k)=0.0;
            else
                mat_A_new(j,k)=A_nume/A_denome;
            end
        end
    end
end
end

% HMM_Viterbi
function vec_hs_seq=HMM_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi)
mat_emission=get_mat_emission(vec_Xi, vec_lambda);

num_of_states=length(mat_A);
num_of_obs=length(vec_Xi);

mat_hs_seq=zeros(num_of_states, num_of_obs);
vec_logp_seq=zeros(num_of_states,1);

% n=1
for j=1:num_of_states
    mat_hs_seq(j,1)=j;
    vec_logp_seq(j)=log(vec_pi(j)*mat_emission(1,j))/log(10);
end

% n>1
for n=2:num_of_obs
    % copy the seq. up to n-1
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
        [max_element,max_pos]=max(vec_h_logprob_i);
        vec_logp_seq_buf(j)=max_element;
        mat_hs_seq_buf(j,:)=mat_hs_seq(max_pos,:);
        mat_hs_seq_buf(j,n)=j;
    end
    % updata the seq.
    mat_hs_seq=mat_hs_seq_buf;
    vec_logp_seq=vec_logp_seq_buf;
end
[max_element, max_pos]=max(vec_logp_seq);
vec_hs_seq=mat_hs_seq(max_pos,:);
end

% get_hmm_ratefunc
function rate_func= get_hmm_ratefunc(spike_time, bin_width)
EMloop_num=5000;		% number of EM itteration
mat_A=[0.999 0.001; 0.001 0.999];
vec_pi=[0.5 0.5];
mean_rate=length(spike_time)/(spike_time(length(spike_time))-spike_time(1));
vec_lambda=[(mean_rate*0.75)*bin_width (mean_rate*1.25)*bin_width];


% 2D array stores
% 1: begining time of each bins in second
% 2: rate of each bin

% get hmm rate func
vec_spkt=zeros(length(spike_time),1);
for i=1:length(spike_time)
    vec_spkt(i)=spike_time(i)-spike_time(1);
end
vec_Xi=get_vec_Xi(vec_spkt, bin_width);

[mat_Gamma, mat_Xi]=HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi);
mat_A_old = mat_A;
vec_pi_old=vec_pi;
vec_lambda_old=vec_lambda;
loop=0;
flag=0;
while (loop<=EMloop_num && flag==0)
    [vec_pi_new, vec_lambda_new, mat_A_new]=HMM_M_step(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi);
    
    vec_pi=vec_pi_new;
    vec_lambda=vec_lambda_new;
    mat_A=mat_A_new;
    
    sum_check=0.0;
    num_state=length(vec_pi);
    
    for i=1:num_state
        for j=1:num_state
            sum_check=sum_check+abs(mat_A_old(i,j)-mat_A(i,j));
        end
        sum_check=sum_check+abs(vec_pi_old(i)-vec_pi(i));
        sum_check=sum_check+abs(vec_lambda_old(i)-vec_lambda(i));
    end
    if sum_check/(1.0*num_state*(num_state+2))<10^(-7)
        flag=flag+1;
    end
    mat_A_old=mat_A;
    vec_pi_old=vec_pi;
    vec_lambda_old=vec_lambda;
    
    [mat_Gamma, mat_Xi]=HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi);
    
    loop=loop+1;
end

vec_hidden=HMM_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi);

rate_func=zeros(length(vec_Xi),2);

%%%%c_time=0.0;
onset = spike_time(1) - 0.001 * (spike_time(length(spike_time)) - spike_time(1));
c_time = onset;
for n=1:length(vec_Xi)
    state_id=vec_hidden(n);
    rate_func(n,1)=round(c_time*100)/100.0;
    rate_func(n,2)=round(vec_lambda(state_id)*100)/(bin_width*100.0);
    c_time=c_time+bin_width;
end
end

function drawHMM(rate_func)
x = rate_func(:, 1);
y = rate_func(:, 2);

ind = 1;
x_new(ind) = x(1);
y_new(ind) = min(y);
ind = ind + 1;
x_new(ind) = x(1);
y_new(ind) = y(1);
ind = ind + 1;

for i = 2 : length(y)
    if y(i - 1) ~= y(i)
        t = (x(i - 1) + x(i)) / 2;
        x_new(ind) = t;
        y_new(ind) = y(i - 1);
        ind = ind + 1;
        x_new(ind) = t;
        y_new(ind) = y(i);
        ind = ind + 1;
    end
end
x_new(ind) = x(length(x));
y_new(ind) = y(length(y));
plot(x_new, y_new);
axis([min(x) max(x) 0 max(y) * 1.1]);
end

