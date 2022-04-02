function SELSFDCorrelatedRealMMSE = Z_Func_LSFD_CorrelatedSMMSE(RPsi, R, DataPowerMatrix,PilotPowerMatrix,L,K,N,tau_p,tau_c,alpha,pilotIndex)

SELSFDCorrelatedRealMMSE = zeros(1,K);
Bmatrix = zeros(L,K,K);
Csmallmatrix =  zeros(L, K, K); 
Dmatrix = zeros(L, K);
PilotPower=PilotPowerMatrix(1,1);

% Compute Bmatrix
for l = 1 : L
    for k = 1:K
        ind=find(pilotIndex==pilotIndex(k,1))';
        for k1=1:length(ind)
            Bmatrix(l,k,ind(k1)) = alpha^3*trace(PilotPower*tau_p*RPsi(:,:,l,k)*R(:,:,l,ind(k1)));
        end
    end
end

% Compute Cmatrix
for l = 1 : L
    for k = 1 : K
        for k1 = 1 : K
            Csmallmatrix(l,k,k1) =  alpha^4*trace(PilotPower*tau_p*R(:,:,l,k1)*R(:,:,l,k)*RPsi(:,:,l,k));
        end
    end
end

% Compute Dmatrix
for l = 1 : L
     temp1=0;
    for k = 1 : K
        for k1=1:K
            temp1=temp1+PilotPower*R(:,:,l,k1);
        end
        Dmatrix(l,k) = trace(alpha^3*(1-alpha)*diag(diag(temp1+eye(N)))*PilotPower*tau_p*RPsi(:,:,l,k)*R(:,:,l,k));
    end
end


ctemp_1 = zeros(L,1);
ctemp_2=0;
% Compute SINR values
for k = 1:K   
    for l = 1: L
        ctemp_1(l,1)=Dmatrix(l,k)+alpha*Bmatrix(l,k,k);
        for k1 = 1: K
            ctemp_1(l,1) =  ctemp_1(l,1)+ DataPowerMatrix(1,k1)*Csmallmatrix(l,k,k1);
        end
    end
    ind=find(pilotIndex==pilotIndex(k,1))';
    for   k2 =1:length(ind)
        if ind(k2)~=k
            btemp_1= squeeze(Bmatrix(:,k,ind(k2)));
            ctemp_2 = ctemp_2+DataPowerMatrix(1,ind(k2))*(btemp_1*btemp_1');
        end
    end
    clk=ctemp_2+diag(ctemp_1);
    btemp= squeeze(Bmatrix(:,k,k));
    SELSFDCorrelatedRealMMSE(k) = (1-tau_p/tau_c)*log2(1+DataPowerMatrix(1,k)*btemp'*(clk\btemp));
end
end
