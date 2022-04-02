function [R,pilotIndex,APpositions,UEpositions,distancesUE] = Z_Z_generateSetup_threeslope(L,K,N,tau_p,nbrOfSetups)
%INPUT:
%L               = Number of APs per setup
%K               = Number of UEs in the network
%N               = Number of antennas per AP
%tau_p           = Number of orthogonal pilots
%nbrOfSetups     = Number of setups with random UE and AP locations
%p               = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%gainOverNoisedB  = Matrix with dimension L x K x nbrOfSetups where
%                   element (l,k,n) is the channel gain (normalized by the
%                   noise variance) between AP l and UE k in setup n, with
%                   correlated shadowing as in [15]
%R                = Matrix with dimension N x N x L x K x nbrOfSetups
%                   where (:,:,l,k,n) is the spatial correlation matrix
%                   between AP l and UE k in setup n, normalized by noise,
%                   with correlated shadowing as in [15]
%pilotIndex     = Matrix with dimension K x nbrOfSetups containing the
%                   pilot assigned to the UEs for cell-free network, using
%                   the greedy pilot assignment in [15] and assuming
%                   correlated shadowing
%APpositions     = Vector of length L with the AP locations, where the real
%                  part is the horizontal position and the imaginary part
%                  is the vertical position
%UEpositions     = Vector of length K with UE positions, measured in the
%                  same way as APpositions
%distancesUE     = Matrix with same dimension as gainOverNoisedB containing
%                  the distances in meter between APs and UEs



%% Define simulation setup

%Size of the coverage area (as a square with wrap-around)
squareLength = 1000; %meter

%Communication bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFigure = 9;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Parameters for the shadow fading from [15]
sigma_sf = 8;
delta = 0.5;
decorr = 100;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Angular standard deviation around the nominal angle (measured in degrees)
ASDdeg = 15;


%Prepare to save results
gainOverNoisedB = zeros(L,K,nbrOfSetups);
R = zeros(N,N,L,K,nbrOfSetups);
distancesUE = zeros(L,K,nbrOfSetups);
pilotIndex = zeros(K,nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Random AP locations with uniform distribution
    APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;
    
    %Random UE locations with uniform distribution
    UEpositions = (rand(K,1) + 1i*rand(K,1)) * squareLength;
    
    
    %Compute alternative AP locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
    UEpositionsWrapped = repmat(UEpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[K 1]);
    
    
    %Compute the correlation matrices for the shadow fading
    shadowCorrMatrix_APs = zeros(L,L);
    shadowCorrMatrix_UEs = zeros(K,K);
    
    for l = 1:L
        distancetoAP = min(abs(APpositionsWrapped - repmat(APpositions(l),size(APpositionsWrapped))),[],2);
        shadowCorrMatrix_APs(:,l) = 2.^(-distancetoAP/decorr);
    end
    
    for k = 1:K
        distancetoUE = min(abs(UEpositionsWrapped - repmat(UEpositions(k),size(UEpositionsWrapped))),[],2);
        shadowCorrMatrix_UEs(:,k) = 2.^(-distancetoUE/decorr);
    end
    
    
    %Generate shadow fading realizations
    a = sigma_sf*sqrtm(shadowCorrMatrix_APs)*randn(L,1);
    b = sigma_sf*sqrtm(shadowCorrMatrix_UEs)*randn(K,1);
    
    for k = 1:K
        
        %Compute distances between each of the APs to UE k
        [distancetoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
        distancesUE(:,k,n) = distancetoUE;
        
        %Compute the channel gain divided by the noise power
        gainOverNoisedB(:,k,n) = pathloss_threeslope(distancesUE(:,k,n)) - noiseVariancedBm;
        
        %Add shadow fading to all channels from APs to UE k that have a
        %distance that is larger than 50 meters
        gainOverNoisedB(distancetoUE>50,k,n) = gainOverNoisedB(distancetoUE>50,k,n) + sqrt(delta)*a(distancetoUE>50) + sqrt(1-delta)*b(k);
        
        
        %Go through all APs
        for l = 1:L
            
            %Compute nominal angle between UE k and AP l
            angletoUE = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l)));
            
            %Generate normalized spatial correlation matrix using the local
            %scattering model
            correlationMatrix = functionRlocalscattering(N,angletoUE,ASDdeg,antennaSpacing);
            R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*correlationMatrix;
            
        end
        
    end
    
    pilotIndex(:,n)=randi(tau_p,K,1);
%     %Assign random pilots while guaranteeing that each pilot is used
%     %equally many times, as an initiation to the greedy algorithm from [15]
%     pilotIndex(:,n) = mod(randperm(K),tau_p)+1;
%     
%     
%     %Extract the normalized large-scale fading coefficients in linear scale
%     betaVal = db2pow(gainOverNoisedB(:,:,n));
    
    
%     %Run the greedy pilot assignment algorithm from [15] for the cell-free
%     %network, with as many iterations as there are UEs
%     for m = 1:K
%         
%         %Compute the SE for all UEs
%         SE_CF = functionComputeSE_CF_uplink_ngo(p_int,p_int*ones(K,1),L,K,tau_p,tau_p+1,pilotIndex(:,n),betaVal);
%         
%         %Find the UE with the lowest SE
%         [~,UEindex] = min(SE_CF);
%         
%         %Remove the pilot assignment for this UE
%         pilotIndex(UEindex,n) = 0;
%         
%         %Compute the interference level for each of the pilots
%         pilotInterference = zeros(tau_p,1);
%         
%         for t = 1:tau_p
%             pilotInterference(t) = sum(sum(betaVal(:,pilotIndex(:,n) == t)));
%         end
%         
%         %Find the pilot with the lowest interference level
%         [~,bestPilot] = min(pilotInterference);
%         
%         %Assign this pilot to the UE
%         pilotIndex(UEindex,n) = bestPilot;
%         
%     end   
end


