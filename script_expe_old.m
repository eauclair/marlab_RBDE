%% Expérimentations sur modèle "simple"
par=12;
n=30;
T=15;
Q=10;
k0=5;
nBlock=5;
%Nombre de paramètres du modèle
nIni=1;%Nombre d'initialisations différentes
nImp=1;%Nombre d'impulseurs différents
nInh=1;%Nombre d'inhibiteurs différents
nAct=1;%Nombre d'actions différentes
nEt=nImp+nInh;%Nombre total d'étiquettes
nSBM=5;
nExp=1;
%Covariable
A=unidrnd(nAct,n,T,Q);
%Types d'initialisation
I=unidrnd(nIni,n,1);
%Blocs aléatoires
SBMblock=unidrnd(nBlock, n, 1);
% Paramètres SBM
XSBM=[0.8 0.3 0.2 0.4 0.4];
%parfor nExp=1:15
    %Paramétrisation 1
    switch mod(nExp,3)
        case 1
            script_param_basique_1;
        case 2
            script_param_basique_2;
        case 0
            script_param_basique_3;
    end
    %Simulation de graphe/données
    script_simu;
    %Gsel = [i j type(inh ou imp) force present/absent]
    Gsel0=zeros(0,5);
    Gsel20=Gr(rand(size(Gr,1),1)<0.2,:);
    dataRBDE={{ini imp inh act} SBMblock Gr A I D Gsel20};
    parsave(strcat('expe_RBDE_simu_data_',num2str(nExp),'.mat'),dataRBDE);
    Gsel=Gsel0;
    script_apprentissage;%Apprentissage de structure
    resAppt_noKn={graphFinTot paramFinTot timeERFinTot lvFinTot};
    parsave(strcat('expe_RBDE_simu_noKn_',num2str(nExp),'.mat'),resAppt_noKn);
    script_apprentissage_SBM;%Apprentissage de structure avec prior SBM
    resAppt_SBM={graphFinTot paramFinTot timeERFinTot lvFinTot paramSBMFinTot};
    parsave(strcat('expe_RBDE_simu_SBM_',num2str(nExp),'.mat'),resAppt_SBM);
    Gsel=Gsel20;
    script_apprentissage;%Apprentissage de structure avec arcs connus
    resAppt_brutKn={graphFinTot paramFinTot timeERFinTot lvFinTot};
    parsave(strcat('expe_RBDE_simu_brutKn_',num2str(nExp),'.mat'),resAppt_brutKn);
%end