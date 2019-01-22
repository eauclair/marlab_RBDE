%function[Ale, Ble, Ri, Gimp, Ginh,Init]=varOptGSBMEsp(D, Gsel, i, k, nIni, nImp, nInh)
function[C, O]=fctSBM(X, XSBM, Dat, A, TL, Ri, Gi, Gfn, nbV, i)
    %n=size(D, 1);
    T=size(Dat, 2);
    nDs=size(Dat,3);
    C=zeros(1,nbV);%Matrice C
    O=0;%R??sultat de la fonction
    e=X(1);
    a=X(2);
    b=X(3);
    c=X(4);
    mu=X(5);
    %Recolonisation (fixe)
    for t=1:T-1
        for s=1:nDs
            O=O+((1-Dat(i,t,s))*Dat(i,t+1,s)*A(t))*log(e);%Recolonisation+Protection
            O=O+((1-Dat(i,t,s))*Dat(i,t+1,s)*(1-A(t)))*log(mu*e);%Recolonisation non protection
            O=O+((1-Dat(i,t,s))*(1-Dat(i,t+1,s))*A(t))*log(1-e);%non recolonisation + protection
            O=O+((1-Dat(i,t,s))*(1-Dat(i,t+1,s))*(1-A(t)))*log(1-mu*e);%non recolonisation non protection
        end
    end
    %Survie (variable)
    for r=1:size(Ri,1)
        %i=Ri(r,1);
        t=Ri(r,1);
        s=Ri(r,7);
        %v=Ri(r,7);
        Pp=1-((1-a)^Ri(r,2));%Quantit?? "1-(1-rhop)^Np
        Pf=1-((1-b)^Ri(r,3));%Quantit?? "1-(1-rhof)^Nf
        Pn=((1-c)^Ri(r,4));%Quantit?? "1-(1-rho-)^N-
        switch Ri(r,5)
            case 0%pf
                P=Pp*Pf*Pn;
            case 1%.f
                P=Pf*Pn;
            case 2%p.
                P=Pp*Pn;
            case 3%..
                P=Pn;
        end
        %P=P+tol;
        if t<T
            if P>0
                if (Dat(i,t,s)*Dat(i,t+1,s)*A(t))%survie, protection
                    C(Ri(r,6))=log(P);
                end
                if (Dat(i,t,s)*Dat(i,t+1,s)*(1-A(t)))%survie, non protection
                   C(Ri(r,6))=log(mu*P);
                end
            end
            if P<1
                if (Dat(i,t,s)*(1-Dat(i,t+1,s))*A(t))%non survie, protection
                    C(Ri(r,6))=log(1-P);
                end
                if (Dat(i,t,s)*(1-Dat(i,t+1,s))*(1-A(t)))%non survie, non protection
                    C(Ri(r,6))=log(1-mu*P);
                end
            end
            %C(r)=q;
        end
    end
          
    if size(XSBM,2)>1
         % la partie de la logvraisemblance qui vient du SBM
        n = size(Dat,1);
        alpha = XSBM(1);
        beta1 = XSBM(2);
        beta2 = XSBM(3);
        alphaft = XSBM(4);
        betaft = XSBM(5);

        for j = 1 :n
            deltaji = TL(j) - TL(i);
            isft=TL(j)*TL(i);
            isft=~isft;
            C(Gi(j,2)) = (alpha*deltaji - log(1-beta2))*(deltaji <0) + (log(alphaft) -log(1-alphaft))*isft ;
            C(Gi(j,3)) = (log(beta2) -log(1-beta2))*(deltaji <=0) +  (log(beta1) - log(1-beta1))*(deltaji>0) + (log(betaft) -log(1-betaft))*isft;
            C(Gfn(j))  = -(log(beta2) -log(1-beta2))*(deltaji <0);

            %ancienne formule, erreur
            %C(Gi(j,2)) = (alpha*deltaji - log(1 + exp(alpha*deltaji)))*(deltaji <0);
            %C(Gi(j,3)) =  (log(beta1) - log(1-beta1)*(deltaji>0)  + (log(beta2) -log(1-beta2)))*(deltaji <=0)  ; 
        end
    end
    
  
    
    %C=C+probsup;
    C=-C;