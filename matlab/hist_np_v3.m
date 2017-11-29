function  OptN  = hist_np_v3( spike_time )
%Function 'hist_np_v3' returns the optimal number of bins of time-histogram. 
%
%Input argument
%spike_time: a vector of spike times {t_i}, spike_time=[t_1,t_2,t_3,...,t_n] (t_1<t_2<t_3< ... <t_n)
%
%Output argument
%optN: the Optimal number of bins.
%
%
%Ref[1]: 
%Takahiro Omi & Shigeru Shinomoto, "Optimizing time histograms for non-Poissonian spike trains", Neural Computation 23, 3125 (2011).
%
%Contact:
%omitakahiro@gmail.com, shinomoto@scphys.kyoto-u.ac.jp
%
%revised by Kazuki Namamura 2017/11/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isrow(spike_time) == 0 
    spike_time = spike_time';
end

min_sp = min(spike_time);    % the time of the first spike
max_sp = max(spike_time);    % the time of the last spike

N = 1:500;                % # of bins vector

D = (max_sp-min_sp)./N;     % bin size  vector
Cost = zeros(1,length(N));  % cost function vector 

SN = 10;                    % # of partitioning positions for shift average

Lv = Calc_Lv( spike_time ); % Lv of the spike train

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computing the cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(N)
    shift = linspace(-D(i)/2,D(i)/2,SN);
    w = zeros(1,SN);
    for p = 1 : SN
        
       edges = linspace(min_sp+shift(p),...
                        max_sp+shift(p),N(i)+1);   % Bin edges

       edges(1)= -inf;
       edges(N(i)+1)= inf;
       
       k=histc(spike_time,edges); 
       k=k(1:N(i));% # of spikes in each bin (Step 1 in ref[1] pp 3129-3130)
       
       %f = arrayfun(@est_fano,k,ones(1,N(i))*2.0*Lv/(3.0-Lv)); % fano factor of each bin (Step 2 in ref[1] pp 3129-3130)
       %2017/11/16 revised.
       f=ones(1,N(i));
       f(k>2)=2.0*Lv/(3.0-Lv);
       
       %Cost(i)= Cost(i) + ( mean(2.0*f.*k)-var(k,1))/ D(i)^2; % cost function for # of bins = N(i) (Step 3,4 in ref[1] pp 3129-3130)
       w(p)=mean(2.0*f.*k)-var(k,1);
       %2017/11/16 revised.
    end
    Cost(i)=sum(w)/D(i)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Minimizing the cost function (Step 5 in ref[1] pp 3129-3130)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,OptN] = min(Cost);
%hist(spike_time,OptN);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function "Calc_Lv" returns the local variation Lv of the spike train
function  Lv  = Calc_Lv( spike_time )
sp1=spike_time(1:length(spike_time)-2);
sp2=spike_time(2:length(spike_time)-1);
sp3=spike_time(3:length(spike_time)  );

Lv = 3.0*mean(((sp3-2.0*sp2+sp1)./(sp3-sp1)).^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate the fano factor for each bin (Step 2 in ref[1] pp %3129-3130)
%this function is not used here.

%function  f = est_fano( count , est_fano )
%if count<3
%    f = 1.0;
%else
%    f = est_fano;
%end
%end

