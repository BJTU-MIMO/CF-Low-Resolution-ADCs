warning('off');
%% Define simulation setup
%Number of Monte Carlo setups
nbrOfSetups = 500;
%Number of APs in the cell-free network
L = 50;
%Number of UEs
K = 30;
%Number of antennas per AP
N = 20;
%Length of the coherence block
tau_c = 200;
%初始值
p_int= 200;
p_int_data=200;
%Number of Iterations for Optimization
NumIter = 200; 
%Number of pilots per coherence block
tau_p = K/5;
%功率
DataPowerMatrix = p_int_data*ones(1,K);
PilotPowerMatrix = p_int*ones(1,K);
% IntPowerMatrix =sqrt(p_int)*ones(1,K);
IntPowerMatrix = sqrt(p_int_data)*rand(1,K);


SumRateLSFDCorrelatedSMMSE= zeros(nbrOfSetups,K);
SumRateOptLSFDCorrelatedSMMSE = zeros(nbrOfSetups,K);

flag=0;

for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    %Generate one setup with UEs at random locations
    [R,pilotIndex] = Z_Z_generateSetup_threeslope(L,K,N,tau_p,1);

    
%     for b=1:1:8
        b=10;
        alpha=finda(b);
        [Rpsi] = Z_Q_functionChannelEstimates(R,L,N,tau_p,pilotIndex,p_int,alpha);
        SELSFDCorrelatedSMMSE = Z_Func_LSFD_CorrelatedSMMSE(Rpsi,R,DataPowerMatrix,PilotPowerMatrix,L,K,N,tau_p,tau_c,alpha,pilotIndex);
        [SEOptLSFDCorrelatedSMMSE, OptDataPowermatrixLFSD_SMMSE,Asave] = Z_Func_OptLSFD_CorrelatedSMMSE(IntPowerMatrix,Rpsi,R,p_int,PilotPowerMatrix,L,K,N,tau_p,tau_c,NumIter,alpha,pilotIndex);
       if p_int*K-sum(OptDataPowermatrixLFSD_SMMSE)<=1
           flag=flag+1;
       end
        SumRateLSFDCorrelatedSMMSE(n,:) = SELSFDCorrelatedSMMSE';
        SumRateOptLSFDCorrelatedSMMSE(n,:) = SEOptLSFDCorrelatedSMMSE;
        clear B;
%     end
end
% Results.MeanSumRateLSFDCorrelatedSMMSE = mean(SumRateLSFDCorrelatedSMMSE);
% Results.MeanSumRateOptLSFDCorrelatedSMMSE =mean(SumRateOptLSFDCorrelatedSMMSE);

figure; 
box on; 
% x=1:1:20;
% With MMSE estimator
p1=cdfplot(sum(SumRateLSFDCorrelatedSMMSE,2));
set(p1,'color','b');
hold on
p2=cdfplot(sum(SumRateOptLSFDCorrelatedSMMSE,2));
set(p2,'color','r');
% xlabel('Number of quantization bits');
xlabel('Sum SE [b/s/Hz]');
hold off
grid on