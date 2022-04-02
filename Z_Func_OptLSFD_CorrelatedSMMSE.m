function [SEOptLSFDCorrelatedSMMSE_r, DataPowermatrix,Asave] =  Z_Func_OptLSFD_CorrelatedSMMSE(IntDataPowerMatrix,RPsi, R, DataPowerMax, PilotPowerMatrix,L,K,N,tau_p,tau_c,NumIter,alpha,pilotIndex)

PilotPower=PilotPowerMatrix(1,1);
SEOptLSFDCorrelatedSMMSE = zeros(NumIter,K);% Space for SINRvalues;
Bmatrix = Compute_Bmatrix(tau_p,PilotPower,RPsi, R,L, K,alpha,pilotIndex);
Csmallmatrix = Compute_Cmatrix(tau_p,PilotPower,RPsi,R,L,K,alpha);
Dmatrix = Compute_Dmatrix(tau_p,PilotPower,R,L,K,N,alpha,RPsi);
SqrtDataPowermatrix = IntDataPowerMatrix;
Amatrix = zeros(L,K);



FlagNote = 0;
for k = 1:K
    Clk  = Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,SqrtDataPowermatrix.^2,L,K,k,alpha, FlagNote,pilotIndex);
    Amatrix(:,k) = squeeze(Clk)\squeeze(Bmatrix(:,k,k));
end
FlagNote = 1;

Asave = zeros(NumIter,1);
DataPowermatrixSave = zeros(K,NumIter);
for Iter = 1:NumIter
    
    % Update U matrix
    [Umatrix, ~, ShareTermForA] = Compute_Umatrix(SqrtDataPowermatrix,Amatrix, Bmatrix,Csmallmatrix,Dmatrix, L, K,pilotIndex,alpha);
    
%     % Update W matrix
%     Wmatrix  = Compute_Wmatrix(SqrtDataPowermatrix,Umatrix,ShareTermForE,Amatrix,Bmatrix, K);
    
    % Update A matrix
    Amatrix  = Compute_Amatrix(SqrtDataPowermatrix,Bmatrix, Csmallmatrix, Dmatrix, ShareTermForA, L, K,alpha,FlagNote,pilotIndex);
    SqrtDataPowermatrix = Compute_Rho(Umatrix,Amatrix,Bmatrix,Csmallmatrix,DataPowerMax,K,L,pilotIndex);
    DataPowermatrix = SqrtDataPowermatrix.^2;
    DataPowermatrixSave(:,Iter) = DataPowermatrix;
    for k = 1:K
        Clk  = Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,DataPowermatrix,L,K,k,alpha, 0,pilotIndex);
        SEOptLSFDCorrelatedSMMSE(Iter,k) = (1-tau_p/tau_c)*log2(1+DataPowermatrix(1,k)*squeeze(Bmatrix(:,k,k))'*((Clk)\squeeze(Bmatrix(:,k,k))));
    end
    Asave(Iter)= sum(SEOptLSFDCorrelatedSMMSE(Iter,:));
end

[~,IterMax] = max(Asave);
SEOptLSFDCorrelatedSMMSE_r= SEOptLSFDCorrelatedSMMSE(IterMax,:);
DataPowermatrix = DataPowermatrixSave(:,IterMax);
end

function Bmatrix = Compute_Bmatrix(tau_p,PilotPower,RPsi,R,L,K,alpha,pilotIndex)

Bmatrix = zeros(L,K,K);
for l = 1 : L
    for k = 1:K
        ind=find(pilotIndex==pilotIndex(k,1))';
        for k1 = 1:length(ind)
            Bmatrix(l,k,ind(k1)) = alpha^3*trace(PilotPower*tau_p*RPsi(:,:,l,k)*R(:,:,l,ind(k1)));
        end
    end
end
end

function Csmallmatrix = Compute_Cmatrix(tau_p,PilotPower,RPsi,R,L,K,alpha)

Csmallmatrix =  zeros(L, K, K); %c_{jk}^it in the paper
% Compute Cmatrix
for l = 1 : L
    for k1 = 1 : K
        for k2 = 1 : K
            Csmallmatrix(l,k1,k2) =  alpha^4*trace(PilotPower*tau_p*R(:,:,l,k2)*R(:,:,l,k1)*RPsi(:,:,l,k1));
        end
    end
end
end

function Dmatrix = Compute_Dmatrix(tau_p,PilotPower,R,L,K,N,alpha,RPsi)

Dmatrix = zeros(L, K);
for l = 1 : L
     temp1=0;
    for k = 1 : K
        for k1=1:K
            temp1=temp1+PilotPower*R(:,:,l,k1);
        end
        Dmatrix(l,k) = trace(alpha^3*(1-alpha)*diag(diag(temp1+eye(N)))*PilotPower*tau_p*RPsi(:,:,l,k)*R(:,:,l,k));
    end
end

end


function Clk  = Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,DataPower,L,K,k,alpha,FlagNote,pilotIndex)

% compute the second part of Clk
ctemp_1 = zeros(L,1);
ctemp_2=0;
% Compute SINR values  
for l = 1: L
    ctemp_1(l,1)=Dmatrix(l,k)+alpha*Bmatrix(l,k,k);
    for k1 = 1: K
        ctemp_1(l,1) =  ctemp_1(l,1)+ DataPower(1,k1)*Csmallmatrix(l,k,k1);
    end
end
ind=find(pilotIndex==pilotIndex(k,1))';
for   k2 =1:length(ind)
    btemp_1= squeeze(Bmatrix(:,k,ind(k2)));
    if ind(k2)==k
        if FlagNote==0
            ctemp_2=0;
        else
            ctemp_2 = ctemp_2+DataPower(1,ind(k2))*(btemp_1*btemp_1');
        end
    else
        ctemp_2 = ctemp_2+DataPower(1,ind(k2))*(btemp_1*btemp_1');
    end
end
Clk=ctemp_2+diag(ctemp_1);
end

function [Umatrix, ShareTermForE, ShareTermForA] = Compute_Umatrix(SqrtDataPowermatrix,Amatrix, Bmatrix,Csmallmatrix,Dmatrix, L, K,pilotIndex,alpha)

Umatrix = zeros(1,K);
ShareTermForE = zeros(1,K);
ShareTermForA = zeros(1,K);
for k = 1: K
    Numerator_ulk = SqrtDataPowermatrix(1,k)*Amatrix(:,k)'*squeeze(Bmatrix(:,k,k));
    Numerator_ulkForA =Amatrix(:,k)'*squeeze(Bmatrix(:,k,k));
    Denominator_ulk = (Amatrix(:,k).^2)'*Dmatrix(:,k)+alpha*(Amatrix(:,k).^2)'*squeeze(Bmatrix(:,k,k));
    for l = 1 : L
        for k1 = 1 :K
            Denominator_ulk = Denominator_ulk +(SqrtDataPowermatrix(1,k1)^2)*(Amatrix(l,k).^2)*Csmallmatrix(l,k,k1);
        end
    end
    ind=find(pilotIndex==pilotIndex(k,1))';
    for   k2 =1:length(ind)
        Denominator_ulk = Denominator_ulk + (SqrtDataPowermatrix(1,ind(k2))^2)*(Amatrix(:,k)'*squeeze(Bmatrix(:,k,ind(k2))))^2;
    end
    Umatrix(1,k) = Numerator_ulk/Denominator_ulk;
    ShareTermForE(1,k) = Denominator_ulk;
    ShareTermForA(1,k) = Numerator_ulkForA/Denominator_ulk;
end
end

% function Wmatrix  = Compute_Wmatrix(SqrtDataPowermatrix,Umatrix,ShareTermForE,Amatrix,Bmatrix, K)
% Wmatrix = zeros(1, K);
% for k = 1 : K
%     elk = (Umatrix(1,k)^2)*ShareTermForE(1,k) - 2*Umatrix(1,k)*SqrtDataPowermatrix(1,k)*abs(Amatrix(:,k)'*Bmatrix(:,k)) + 1;
%     Wmatrix(1,k) = 1/elk;
% end
% end

function  Amatrix  = Compute_Amatrix(SqrtDataPowermatrix,Bmatrix, Csmallmatrix, Dmatrix, ShareTermForA, L, K,alpha,FlagNote,pilotIndex)
Amatrix = zeros(L,K);
for k = 1 : K
    Clk =  Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,SqrtDataPowermatrix.^2,L,K,k,alpha,FlagNote,pilotIndex);
    Amatrix(:,k) =(1/ShareTermForA(1,k))*(Clk\squeeze(Bmatrix(:,k,k)));
end
end

function SqrtDataPowermatrix = Compute_Rho(Umatrix,Amatrix,Bmatrix,Csmallmatrix,DataPowerMax,K,L,pilotIndex)
SqrtDataPowermatrix = zeros(1,K);
for k = 1 : K
    Numerator_rholk = Umatrix(1,k)*abs(Amatrix(:,k)'*squeeze(Bmatrix(:,k,k)));
    Denominator_rholk = 0;
%     for l=1:L
%         Denominator_rholk = Denominator_rholk + (Umatrix(1,k)^2)*(Amatrix(l,k).^2)*sum(Csmallmatrix(l,k,:),3) ; 
%     end
%     ind=find(pilotIndex==pilotIndex(k,1))';
%     for   k1 =1:length(ind)
%         Denominator_rholk = Denominator_rholk + (Umatrix(1,k)^2)*(Amatrix(:,k)'*squeeze(Bmatrix(:,k,ind(k1))))^2; 
%     end
    for l=1:L
        Denominator_rholk = Denominator_rholk + (Umatrix(1,k)^2)*(Amatrix(l,k).^2)*Csmallmatrix(l,k,k) ; 
    end
    Denominator_rholk = Denominator_rholk + (Umatrix(1,k)^2)*(Amatrix(:,k)'*squeeze(Bmatrix(:,k,k)))^2;
    SqrtDataPowermatrix(1,k) = min(Numerator_rholk/Denominator_rholk,sqrt(DataPowerMax));
end
end
