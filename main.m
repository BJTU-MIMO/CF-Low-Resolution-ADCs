%% Define simulation setup
%Number of Monte Carlo setups
%蒙特卡洛仿真的次数
nbrOfSetups = 300;

%Number of channel realizations per setup
%每次仿真信道实现的次数
nbrOfRealizations = 1000;

%Number of APs in the cell-free network
%AP的数目
L = 100;

%Number of UEs
%用户的数目
K = 40;

%Number of antennas per AP
%每个AP上天线的数目
N = 4;

%Length of the coherence block
%一个相干资源块的长度
tau_c = 200;

%Number of pilots per coherence block
%导频的个数
tau_p = K/4;

%Uplink transmit power per UE (mW)
%发送功率，上行
p = 100;

SE_MR_L2_Q_dis=zeros(K,10,nbrOfSetups);
SE_MMSE_L2_Q_dis=zeros(K,10,nbrOfSetups);

SE_MR_L3_Q_dis=zeros(K,10,nbrOfSetups);
SE_MMSE_L3_Q_dis=zeros(K,10,nbrOfSetups);

SE_MR_L1_Q_dis=zeros(K,10,nbrOfSetups);
SE_MMSE_L1_Q_dis=zeros(K,10,nbrOfSetups);

SE_MR_L2_dis=zeros(K,nbrOfSetups);
SE_MMSE_L2_dis=zeros(K,nbrOfSetups);

SE_MR_L3_dis=zeros(K,nbrOfSetups);
SE_MMSE_L3_dis=zeros(K,nbrOfSetups);

SE_MR_L1_dis=zeros(K,nbrOfSetups);
SE_MMSE_L1_dis=zeros(K,nbrOfSetups);
%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs at random locations
    %随机撒点生成用户和AP的位置，输出信道的相关矩阵R和导频分配向量pilotIndex
    [gainOverNoisedB,R,pilotIndex] = Z_generateSetup_threeslope(L,K,N,tau_p,1,p);
    
    for b=1:10
        %根据量化比特数目确定alpha
        alpha=finda(b);
        %Generate channel realizations, channel estimates, and estimation
        %error correlation matrices for all UEs to the APs 
        %输出信道估计
        [Hhat,H,B] = Z_Q_functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p,alpha);
        %Compute SE using Monte Carlo simulations
        %计算SE
        [SE_MR_Q,SE_MMSE_Q] = Z_Q_ceshi1_functionComputeSE_AP_uplink(Hhat,H,R,B,tau_c,tau_p,nbrOfRealizations,N,K,L,p,alpha);
        SE_MR_L1_Q_dis(:,b,n)=log2(1+SE_MR_Q(:,1));
        SE_MMSE_L1_Q_dis(:,b,n)=log2(1+SE_MMSE_Q(:,1));
        SE_MR_L2_Q_dis(:,b,n)=log2(1+SE_MR_Q(:,2));
        SE_MMSE_L2_Q_dis(:,b,n)=log2(1+SE_MMSE_Q(:,2));
        SE_MR_L3_Q_dis(:,b,n)=log2(1+SE_MR_Q(:,3));
        SE_MMSE_L3_Q_dis(:,b,n)=log2(1+SE_MMSE_Q(:,3));

        %Remove large matrices at the end of analyzing this setup
        clear B H Hhat;
    end
    %idealADC情况
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the APs 
    [Hhat,H,B] = Z_functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    %Compute SE using Monte Carlo simulations
    [SE_MR,SE_MMSE] = Z_functionComputeSE1_AP_uplink(Hhat,H,R,B,tau_c,tau_p,nbrOfRealizations,N,K,L,p);
    clear B H Hhat;
    SE_MR_L2_dis(:,n)=log2(1+SE_MR(:,2));
    SE_MMSE_L2_dis(:,n)=log2(1+SE_MMSE(:,2));
    SE_MR_L3_dis(:,n)=log2(1+SE_MR(:,3));
    SE_MMSE_L3_dis(:,n)=log2(1+SE_MMSE(:,3));
    SE_MR_L1_dis(:,n)=log2(1+SE_MR(:,1));
    SE_MMSE_L1_dis(:,n)=log2(1+SE_MMSE(:,1));
end

r_MR_Q=zeros(3,10);
r_MMSE_Q=zeros(3,10);

%求均值
for x=1:10
    r_MR_Q(2,x)=mean(mean(squeeze(SE_MR_L2_Q_dis(:,x,:))));
    r_MMSE_Q(2,x)=mean(mean(squeeze(SE_MMSE_L2_Q_dis(:,x,:))));
    r_MR_Q(3,x)=mean(mean(squeeze(SE_MR_L3_Q_dis(:,x,:))));
    r_MMSE_Q(3,x)=mean(mean(squeeze(SE_MMSE_L3_Q_dis(:,x,:))));
    r_MR_Q(1,x)=mean(mean(squeeze(SE_MR_L1_Q_dis(:,x,:))));
    r_MMSE_Q(1,x)=mean(mean(squeeze(SE_MMSE_L1_Q_dis(:,x,:))));
    
end


r_MR_L2=mean(mean(SE_MR_L2_dis))*ones(1,10);
r_MMSE_L2=mean(mean(SE_MMSE_L2_dis))*ones(1,10);
r_MR_L3=mean(mean(SE_MR_L3_dis))*ones(1,10);
r_MMSE_L3=mean(mean(SE_MMSE_L3_dis))*ones(1,10);
r_MR_L1=mean(mean(SE_MR_L1_dis))*ones(1,10);
r_MMSE_L1=mean(mean(SE_MMSE_L1_dis))*ones(1,10);


figure(1);
hold on; box on;
x=1:1:10;
plot(x,r_MMSE_Q(3,:),'b-','LineWidth',2);
plot(x,r_MMSE_L3,'b--','LineWidth',2);
plot(x,r_MR_Q(3,:),'r-','LineWidth',2);
plot(x,r_MR_L3,'r--','LineWidth',2);

plot(x,r_MMSE_Q(2,:),'g-','LineWidth',2);
plot(x,r_MMSE_L2,'g--','LineWidth',2);
plot(x,r_MR_Q(2,:),'k-','LineWidth',2);
plot(x,r_MR_L2,'k--','LineWidth',2);

plot(x,r_MMSE_Q(1,:),'m-','LineWidth',2);
plot(x,r_MMSE_L1,'m--','LineWidth',2);
plot(x,r_MR_Q(1,:),'c-','LineWidth',2);
plot(x,r_MR_L1,'c--','LineWidth',2);


ylabel('Average per-user spectral efficiency ');
xlabel('Quantization bits');
legend({'L3Q (L-MMSE)','L3 (L-MMSE)','L3Q (MR)','L3 (MR)','L2Q (L-MMSE)','L2 (L-MMSE)','L2Q (MR)','L2 (MR)','L1Q (L-MMSE)','L1 (L-MMSE)','L1Q (MR)','L1 (MR)'},'Location','NorthWest');


