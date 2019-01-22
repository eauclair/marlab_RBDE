%% Apprentissage R-E avec prior SBM pour réseau écologique
%Demande : k0 : degré max ;   D : Données présence/absence ; 
%I : Initialisation ; par0, parEq, parKn : paramètres ;
% SBMBlock : bloc du SBM
n=size(D,1);
T=size(D,2);
Q=size(D,3);
%I : Initialisation ; 
%Graphe vide
G=zeros(0,4);
[iniE,impE,inhE,actE]=parEstKn(D,G,A,I,par0,parKn,parEq);
lv0=lvTot(D,G,A,I,iniE,impE,inhE,actE);
%Estimation des paramètres SBM
XSBM0 = [0.01 0.01 0.01 0.01 0.01];
%Changement en matrices d'adjacence
AdjImp=ltoA(G(G(:,3)==1,[1 2 4]),n,1);
AdjInh=ltoA(G(G(:,3)==2,[1 2 4]),n,1);
%estimation de alpha
fsbm = @(alpha)-lvSBMalpha(AdjImp, SBMblock, alpha);
alpha = fmincon(fsbm,XSBM0(1), [],[],[],[],0,1,[],optimset('Display','off')); 
%Estimation de Beta1
beta1 = estimbeta1(AdjInh, SBMblock);
if beta1 == 0 beta1 = 0.1; end
%Estimation de Beta2
beta2 = estimbeta2(AdjImp,AdjInh, SBMblock);
if beta2 == 0 beta2 = 0.1; end
%Estimation de alpha et beta fourre-tout
[alphaft, betaft]=estimft(AdjImp, AdjInh, SBMblock);
if alphaft == 0 alphaft = 0.1; end
if betaft == 0 betaft = 0.1; end
%Calcul de la vraisemblance
XSBM_E = [alpha, beta1, beta2, alphaft, betaft];
lv0 = lv0 +  lvSBM1(AdjImp, AdjInh, SBMblock, XSBM_E);
%% Apprentissage
%Ajout de Cplex
addpath(genpath('/opt/ibm/ILOG/CPLEX_Studio1261/cplex/matlab'));
%Lancement de la parallelisation (si pas encore fait)
if isempty(gcp('nocreate'))
    parpool(par);
end
%Gsel = [i j type(inh ou imp) force present/absent]
Gsel=zeros(0,5);
go=1;
timeER=tic;
while go
    %ILP
    lvC=0;
    GappCell=cell(n,1);
    IappCell=cell(n,1);
    %Contraintes de base pour décrire les 
    [Ale, Ble, indVar0, Ri0]=varILP_core(D, k0, nIni, nImp, nInh);
    parfor i=1:n
        GselI=Gsel(Gsel(:,2)==i,:);
        %Variables de l'ilp (structure de base)
        k=k0;
        Ai=Ale;
        Bi=Ble;
        indVar=indVar0;
        Ri=Ri0;
        if size(GselI,1)
            [A1, B1]=varILP_Gsel(Ai, indVar, GselI);
            Ai=[Ai; A1];Bi=[Bi B1];
        end
        %Variables de l'ilp (contraintes structurelles)
        [A1, B1]=varILP_structC(Ai, indVar,i, nImp, nInh);
        Ai=[Ai; A1];Bi=[Bi B1];
        [A1, B1]=varILP_foodWeb(Ai, indVar);
        Ai=[Ai; A1];Bi=[Bi B1];
        nbV=size(Ai,2);
        %Fonction objectif ILP
        %Ajout des bornes inf et sup comme contraintes
        Ab=speye(nbV);Bb=ones(nbV,1);
        Ab=[Ab;-speye(nbV)];Bb=[Bb;zeros(nbV,1)];
        Ao=[Ai;Ab];Bo=[Bi';Bb];
        Co=lvILP(D, A, i, iniE, impE, inhE, actE, nbV, Ri);
        Co=lvILP_SBM1(D, i, Co, SBMblock', XSBM_E, indVar);
        [GappI, IappI]=restoration_step_ILP(Ao, Bo, -Co, indVar, i);
        GappCell{i}=GappI;
        IappCell{i}=IappI;
    end
    Gapp=cell2mat(GappCell);
    Iapp=cell2mat(IappCell);
    %Estimation des paramètres
    [iniE,impE,inhE,actE]=parEstKn(D,Gapp,A,Iapp,par0,parKn,parEq);
    %Calcul de la vraisemblance
    lv1=lvTot(D,Gapp,A,Iapp,iniE,impE,inhE,actE);
    %Log-vraisemblance : partie SBM
    AdjImp=ltoA(Gapp(find(Gapp(:,3)==1),[1 2 3]),n,1);%penser à changer Gr en G
    AdjInh=ltoA(Gapp(find(Gapp(:,3)==2),[1 2 4]),n,1);
    fsbm = @(alpha)-lvSBMalpha(AdjImp, SBMblock, alpha);
    alpha = fmincon(fsbm,XSBM0(1), [],[],[],[],0,1,[],optimset('Display','off')); 

    beta1 = estimbeta1(AdjInh, SBMblock);
    if beta1 == 0 beta1 = 0.1; end

    beta2 = estimbeta2(AdjImp,AdjInh, SBMblock);
    if beta2 == 0 beta2 = 0.1; end

    [alphaft, betaft]=estimft(AdjImp, AdjInh, SBMblock);
    if alphaft == 0 alphaft = 0.1; end
    if betaft == 0 betaft = 0.1; end

    XSBM_E = [alpha, beta1, beta2, alphaft, betaft];
    lv1 = lv1 +  lvSBM1(AdjImp, AdjInh, SBMblock, XSBM_E);
    if lv1>lv0
        G=Gapp;
        I=Iapp;
        lv0=lv1;
    else
        go=0;
    end

end
timeERFinTot=toc(timeER);
graphFinTot=G;%Graphe final
IniFinTot=I;%Initialisation finale
lvFinTot=lv0;%Vraisemblance finale
paramFinTot=[iniE,impE,inhE,actE];%Paramètres finaux
paramSBMFinTot=XSBM_E;%Paramètres finaux SBM