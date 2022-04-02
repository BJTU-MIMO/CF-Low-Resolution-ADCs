function [SE_MR_Q,SE_MMSE_Q] = Z_Q_ceshi1_functionComputeSE_AP_uplink(Hhat,H,R,B,tau_c,tau_p,nbrOfRealizations,N,K,L,p,alpha)
%INPUT:
%Hhat              = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel from
%                    all BSs to UE k at channel realization n.
%H                 = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the true collective channel from all
%                    BSs to UE k at channel realization n.
%R                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix between BS
%                    l and UE k in setup n, normalized by the noise power
%B                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix of the
%                    estimate between BS l and UE k in setup n, normalized
%                    by the noise power
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences and number of UEs per cell
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Total number of UEs
%L                 = Number of APs
%p                 = Matrix K x 1 where element k is the uplink transmit
%                    power of UE k (If it is a scalar, the same value is
%                    used for all users)
%p1                = (Optional) Same as p but only for Level 1
%
%OUTPUT:
%SE_MR     = K x 4 matrix where the (k,n):th element is the uplink SE of 
%            UE k achieved with MR combining at cooperation level n
%SE_MMMSE  = Same as SE_MR but with MMSE or L-MMSE combining

%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
   p = p*ones(K,1);
end

%Store identity matrices of different sizes
eyeN = eye(N);
eyeLN = eye(L*N);

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results,level2&level3&levelz
SE_MR_Q = zeros(K,3);
SE_MMSE_Q = zeros(K,3);

%估计误差
C_tot=zeros(N,N,L);
for k = 1:K
    C_tot = C_tot + (R(:,:,:,k)-B(:,:,:,k));
end
%量化误差，diag
Q_tot=zeros(N,N,L);
for l = 1:L
    for k=1:K
         Q_tot(:,:,l) = Q_tot(:,:,l) + p(k)*R(:,:,l,k);
    end
    Q_tot(:,:,l)= alpha*(1-alpha)*diag(diag((Q_tot(:,:,l)+eyeN)));
end

Nn=sqrt(0.5)*(randn(N,nbrOfRealizations,L) + 1i*randn(N,nbrOfRealizations,L));

%Diagonal matrix with transmit powers and its square root
Dp = diag(p*alpha^2);
Dpsq = diag(sqrt(p)*alpha);

%Prepare to save simulation results
SE_MR_level1 = zeros(K,L);
SE_MMSE_level1 = zeros(K,L);

signal_MR_level23 = zeros(L,K);
scaling_MR_level23 = zeros(L,K);
Gp1_MR_level23 = zeros(L,L,K);
Gp2_MR_level23 = zeros(L,L,K);

signal_MMSE_level23 = zeros(L,K);
scaling_MMSE_level23 = zeros(L,K);
G1_MMSE_level23 = zeros(L,L,K);
G2_MMSE_level23 = zeros(L,L,K);

%% Go through all channel realizations
for n = 1:nbrOfRealizations
  
   %Level 2&3  UatF
    gp1_MR_level2 = zeros(L,K,K);
    gp2_MR_level2 = zeros(L,K);
    gp1_MMSE_level2 = zeros(L,K,K);
    gp2_MMSE_level2 = zeros(L,K);

    %Go through all APs
    for l = 1:L
        
        %level 2&3
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(l-1)*N:l*N,n,:),[N K]);
        
        %计算量化噪声
        Nq = sqrtm(Q_tot(:,:,l))*squeeze(Nn(:,n,l));
        %Compute combining
        V_MR = Hhatallj;
        
        %Go through all UEs
        for k = 1:K
            %%MR combining
            v = V_MR(:,k); %Extract combining vector
            %Level 2&3 
            signal_MR_level23(l,k) = signal_MR_level23(l,k) + alpha*(v'*Hallj(:,k))/nbrOfRealizations;
            gp1_MR_level2(l,:,k) = gp1_MR_level2(l,:,k) + (v'*Hallj)*Dpsq;
            gp2_MR_level2(l,k) = gp2_MR_level2(l,k) + v'*Nq;
            scaling_MR_level23(l,k) = scaling_MR_level23(l,k) + alpha^2 * norm(v).^2/nbrOfRealizations;
            %level 1
            %Compute numerator and denominator of instantaneous SINR
            numerator =p(K)*alpha^2*abs(v'*Hhatallj(:,k))^2;
            denominator = norm((v'*Hhatallj)*Dpsq)^2 + v'*(p(K)*alpha^2*C_tot(:,:,l)+alpha^2*eyeN++Q_tot(:,:,l))*v - numerator;
            %Compute instantaneous SE for one channel realization
            SE_MR_level1(k,l) = SE_MR_level1(k,l) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
            
        
            %%MMSE combining
            v =((Hhatallj*Dp*Hhatallj')+Q_tot(:,:,l)+p(K)*alpha^2*C_tot(:,:,l)+alpha^2*eyeN)\(Hhatallj(:,k)*p(k)*alpha^2);%Extract combining vector
            %Level 2&3 
            signal_MMSE_level23(l,k) = signal_MMSE_level23(l,k) + alpha*(v'*Hallj(:,k))/nbrOfRealizations;
            gp1_MMSE_level2(l,:,k) = gp1_MMSE_level2(l,:,k) + (v'*Hallj)*Dpsq;
            gp2_MMSE_level2(l,k) = gp2_MMSE_level2(l,k) + v'*Nq;
            scaling_MMSE_level23(l,k) = scaling_MMSE_level23(l,k) + alpha^2 * norm(v).^2/nbrOfRealizations;
            %Level 1
            %Compute numerator and denominator of instantaneous SINR
            numerator =p(K)*alpha^2*abs(v'*Hhatallj(:,k))^2;
            denominator = norm((v'*Hhatallj)*Dpsq)^2 + v'*(p(K)*alpha^2*C_tot(:,:,l)+alpha^2*eyeN++Q_tot(:,:,l))*v - numerator;
            %Compute instantaneous SE for one channel realization
            SE_MMSE_level1(k,l) = SE_MMSE_level1(k,l) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        end  
    end
    
    %Compute averaging of terms for Level 2&3
    for k = 1:K
        Gp1_MR_level23(:,:,k) = Gp1_MR_level23(:,:,k) + gp1_MR_level2(:,:,k)*gp1_MR_level2(:,:,k)'/nbrOfRealizations;
        G1_MMSE_level23(:,:,k) = G1_MMSE_level23(:,:,k) + gp1_MMSE_level2(:,:,k)*gp1_MMSE_level2(:,:,k)'/nbrOfRealizations;
        Gp2_MR_level23(:,:,k) = Gp2_MR_level23(:,:,k)+gp2_MR_level2(:,k)*gp2_MR_level2(:,k)'/nbrOfRealizations;
        G2_MMSE_level23(:,:,k) = G2_MMSE_level23(:,:,k)+gp2_MMSE_level2(:,k)*gp2_MMSE_level2(:,k)'/nbrOfRealizations;
    end
    
end

%Compute SE for Level 1
SE_MR_Q(:,1) = max(SE_MR_level1,[],2);
SE_MMSE_Q(:,1) = max(SE_MMSE_level1,[],2);
%Compute SE for Level 2&3
for k = 1:K   
    %With MR combining
    b = signal_MR_level23(:,k);
    A = Gp1_MR_level23(:,:,k) + Gp2_MR_level23(:,:,k) + diag(scaling_MR_level23(:,k)) - p(k)*alpha*(b*b');   
    SE_MR_Q(k,2) = prelogFactor*real(log2(1+p(k)*abs(mean(b)).^2 / mean(mean(A))));  
    SE_MR_Q(k,3) = prelogFactor*real(log2(1+p(k)*b'*(A\b)));  

    %With L-MMSE combining
    b = signal_MMSE_level23(:,k);
    A = G1_MMSE_level23(:,:,k)  + G2_MMSE_level23(:,:,k)+ diag(scaling_MMSE_level23(:,k)) - p(k)*alpha*(b*b');
    SE_MMSE_Q(k,2) = prelogFactor*real(log2(1+p(k)*abs(mean(b)).^2 / mean(mean(A))));  
    SE_MMSE_Q(k,3) = prelogFactor*real(log2(1+p(k)*b'*(A\b)));  
    
end
end



