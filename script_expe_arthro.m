%% Expérimentations sur données "arthropodes"
load('Darth.mat')
load('size.mat');
load('Di.mat')
%Données 'intérieures' (points 2 et 32)
Dint(:,1,:)=Dpa(:,2,:)&Dpa(:,3,:);Dint(:,2,:)=Dpa(:,5,:)&Dpa(:,6,:);
DBi=Dint(:,:,find(lCult==1));%Betterave
DMi=Dint(:,:,find(lCult==2));%Maïs
DSRi=Dint(:,:,find(lCult==3));%Colza été
DWRi=Dint(:,:,find(lCult==4));%Colza hiver
par=12;
n=size(D,1);
T=size(D,2);
Q=size(D,3);
k0=5;
nBlock=4;
%Nombre de paramètres du modèle
nIni=2;%Nombre d'initialisations différentes
nImp=1;%Nombre d'impulseurs différents
nInh=1;%Nombre d'inhibiteurs différents
nAct=0;%Nombre d'actions différentes
nEt=nImp+nInh;%Nombre total d'étiquettes
nSBM=5;
%Covariable
A=unidrnd(nAct,n,T,Q);
%Types d'initialisation
I=unidrnd(nIni,n,1);
%Blocs aléatoires
SBMblock=unidrnd(nBlock, n, 1);
% Paramètres SBM
XSBM=[0.8 0.3 0.2 0.4 0.4];
%Ajout de Cplex
addpath(genpath('/opt/ibm/ILOG/CPLEX_Studio1261/cplex/matlab'));
%Lancement de la parallelisation (si pas encore fait)
if isempty(gcp('nocreate'))
    parpool(par);
end
parfor nExp=1:15
    %Paramétrisation 1
    switch mod(nExp,3)
        case 1
            parExist=[1 0 0 0;0 1 1 0];
        case 2
            parExist=[0 1 1 0;1 0 0 0];
        case 0
            parExist=[1 1 1 0;1 1 1 0];
    end
    [ini, imp, inh, act, par0, parKn, parEq]=dataGen_param(nIni, nImp, nInh, nAct, parExist);
    parVal={ini imp inh act};
%Simulation de graphe/données
    %Graphe simulé
    Gr=GSimulSBM(XSBM,SBMblock',k0);
    %Valeurs d'initialisation (aléatoires)
    %---Etat initial
    t0=zeros(n,Q);
    for i=1:n
        for s=1:Q
            if ini(1,I(i))
                t0(i,s)=rand()<ini(1,I(i));
            else
                t0(i,s)=rand()<ini(2,I(i));
            end
        end
    end
    %Simulation de données
    D=dSimul(Gr,A,t0,I,ini,imp,inh,act);
    %Calcul de log vraisemblance "réelle"
    lvR=lvTot(D,Gr,A,I,ini,imp,inh,act);
    %Gsel = [i j type(inh ou imp) force present/absent]
    Gsel0=zeros(0,5);
    Gsel20=Gr(rand(size(Gr,1),1)<0.2,:);
    Gsel20(:,5)=1;
    dataRBDE={parVal SBMblock Gr A I D Gsel20};
 % Apprentissage
    parsave(strcat('expe_RBDE_simu_data_',num2str(nExp),'.mat'),dataRBDE);
    Gsel=Gsel0;
    [graphFinTot, ~, timeERFinTot, lvFinTot, paramFinTot]=ERlearning(D,A,I,k0,par0,parKn,parEq,Gsel);
    resAppt_noKn={graphFinTot paramFinTot timeERFinTot lvFinTot};
    parsave(strcat('expe_RBDE_simu_noKn_',num2str(nExp),'.mat'),resAppt_noKn);
    [graphFinTot, IniFinTot, timeERFinTot, lvFinTot, paramFinTot, paramSBMFinTot]=ERlearningSBM(D,A,I,k0,par0,parKn,parEq,Gsel,SBMblock);
    resAppt_SBM={graphFinTot paramFinTot timeERFinTot lvFinTot paramSBMFinTot};
    parsave(strcat('expe_RBDE_simu_SBM_',num2str(nExp),'.mat'),resAppt_SBM);
    Gsel=Gsel20;
    [graphFinTot, ~, timeERFinTot, lvFinTot, paramFinTot]=ERlearning(D,A,I,k0,par0,parKn,parEq,Gsel);
    resAppt_brutKn={graphFinTot paramFinTot timeERFinTot lvFinTot};
    parsave(strcat('expe_RBDE_simu_brutKn_',num2str(nExp),'.mat'),resAppt_brutKn);
end