%lvILP(D, A,i,ini, imp, inh, act, nbV)
function[C]=lvILP(D, A,i,ini, imp, inh, act, nbV, Ri)
    C=zeros(1,nbV);%Matrice C
    for r=1:size(Ri,1)
        t=Ri(r,2);
        s=Ri(r,3);
        AppImp=1;
        SurImp=1;
        AppInh=1;
        SurInh=1;
        ind=4;
        %Proba d'échec de tous les impulseurs
        for rimp=1:size(imp,2)
            AppImp=AppImp*(1-imp(1,rimp))^Ri(r,ind+rimp);
            SurImp=SurImp*(1-imp(2,rimp))^Ri(r,ind+rimp);
        end
        ind=ind+size(imp,2);
        %Proba d'échec de tous les inhibiteurs
        for rinh=1:size(inh,2)
            AppInh=AppInh*(1-inh(1,rinh))^Ri(r,ind+rinh);
            SurInh=SurInh*(1-inh(2,rinh))^Ri(r,ind+rinh);
        end
        %Proba d'initialisation (spontanée)
        AppIni=ini(1,Ri(r,4));
        SurIni=ini(2,Ri(r,4));
        %Pénalisation due à l'action
        AppAct=act(1,A(i,t,s));
        SurAct=act(2,A(i,t,s));
        %P=P+tol;
        Papp=AppAct*(AppIni+(1-AppIni)*(1-AppImp))*AppInh;
        Psur=SurAct*(SurIni+(1-SurIni)*(1-SurImp))*SurInh;
        P=(1-D(i,t,s))*Papp+D(i,t,s)*Psur;
        if D(i,t+1,s)
            C(Ri(r,1))=log(P);
        else
            C(Ri(r,1))=log(1-P);
        end
    end
    %C=C+probsup;
    %C=-C;