function [Hhat,H,B] = Z_Q_functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p,alpha)
%INPUT:
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k 
%                    in setup n, normalized by the noise power
%nbrOfRealizations = Number of channel realizations
%L                 = Number of APs
%K                 = Number of UEs in the network
%N                 = Number of antennas per AP
%tau_p             = Number of orthogonal pilots
%pilotIndex        = Vector containing the pilot assigned to each UE
%p                 = Uplink transmit power per UE (same for everyone)
%alpha             = 量化的比特数
%Rq                = 量化噪声的方差
%
%OUTPUT:
%Hhat         = Matrix with dimension L*N x nbrOfRealizations x K where
%               (:,n,k) is the estimated collective channel to UE k at
%               channel realization n.
%H            = Matrix with dimension L*N x nbrOfRealizations x K with the
%               true channel realizations. The matrix is organized in the
%               same way as Hhat_MMSE.
%B            = Matrix with dimension N x N x L x K where (:,:,l,j) is the
%               spatial correlation matrix of the estimate between AP l and
%               UE k in setup n, normalized by the noise power


%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K));


%Go through all channels and apply the spatial correlation matrices
for l = 1:L
    
    for k = 1:K
        
        %Apply correlation to the uncorrelated channel realizations
        Rsqrt = sqrtm(R(:,:,l,k));
        H((l-1)*N+1:l*N,:,k) = sqrt(0.5)*Rsqrt*H((l-1)*N+1:l*N,:,k);
        
    end
    
end


%% Perform channel estimation

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));


%Prepare to store results
Hhat = zeros(L*N,nbrOfRealizations,K);

if nargout>2
    B = zeros(size(R));
end

Rq = zeros(size(R));

%Go through all APs
for l = 1:L
    
    %Go through all pilots
    for t = 1:tau_p
        
        %Compute processed pilot signal for all UEs that use pilot t
        yp = sqrt(p*tau_p)*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + Np(:,:,l,t);
        
        %yp的期望
        Eyp =(p*tau_p*sum(R(:,:,l,t==pilotIndex),4) + eyeN);
        
        %量化噪声的方差
        Rq(:,:,l,t)=alpha*(1-alpha)*diag(diag(Eyp));
        Npq=sqrtm(squeeze(Rq(:,:,l,t)))*Np(:,:,l,t);
        
        %量化
        ypq=alpha*yp+Npq;
        
        %需要求逆的地方
        PInv=alpha^2*Eyp+Rq(:,:,l,t);
        
        %Go through all UEs that use pilot t
        for k = find(t==pilotIndex)'
            
            %Compute the MMSE estimate
            RPsi =R(:,:,l,k) / PInv;
            Hhat((l-1)*N+1:l*N,:,k) =  alpha*sqrt(p*tau_p)*RPsi*ypq;
            
            %Compute the spatial correlation matrix of the estimate
            if nargout>2
                B(:,:,l,k) = alpha^2*p*tau_p*RPsi*R(:,:,l,k);
            end
            
        end
        
    end
    
end


