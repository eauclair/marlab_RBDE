%% Script d'exemple d'utilisation de RBDE
%% Variables communes
%-----------------------------------------------------
%cd '/home/eauclair/MIAT/Work/Code/code_LDBN_FoodWeb'
clear()
par=8;
n=5;
T=15;
Q=5;
k0=4;
%Statut de protection (ici, aucun)
A=ones(n,T,Q);
%Block équilibrés (pour SBM) avec 3 FT
block=[1 1 2 2 2]';
%Block=[zeros(6,1);ones(6,1);ones(6,1)+1;ones(6,1)+2;ones(6,1)+4];
%Block=[zeros(3,1);ones(3,1);ones(3,1)+1;ones(3,1)+2;ones(3,1)+4];

nIni=2;%Nombre d'initialisations différentes
nImp=1;%Nombre d'impulseurs différents
nInh=1;%Nombre d'inhibiteurs différents
nAct=1;%Nombre d'actions différentes
nEt=nImp+nInh;%Nombre total d'étiquettes
nSBM=5;
%Calcul du nombre de variables
nVarR=0;
for h=0:k0
    nVarR=nVarR+(factorial(nEt+h-1)/(factorial(h)*factorial(nEt-1)));
end
nVar=n*nEt+nIni+Q*(T+1)*(3*nEt*(k0+1))+Q*T*nVarR*nIni;
%save('example_varComm.mat', 'A', 'Block', 'k0', 'nIni', 'nImp', 'nInh', 'nAct')
%% Paramètres fixés
%---Type d'initialisation de chaque espèce
I_F=[zeros(8,1);ones(7,1)]+1;
%---Parametres du DBN (utilises pour simuler les donnees)
ini_F=[0.8 0.8;1 0];%
imp_F=[0 0;0.8 0.8];
inh_F=[0;0.8];
act_F=[1;1];
%--Paramètres SBM fixés
XSBM_F=[0.8 0.5 0.3 0.4 0.4];
%% Paramètres générés aléatoirement
%---Type d'initialisation de chaque espèce
%I_R=unidrnd(nIni,n,1);
%---Param??tres du RBDE
ini_R=rand(2,nIni);
imp_R=rand(2,nImp);
inh_R=rand(2,nInh);
act_R=[1;1];
%--Paramètres du SBM
XSBM_R=rand(1,nSBM);%XSBM : [alpha beta1 beta2 alphaft betaft]
%% Pas d'inhibiteurs
XSBM_R([2 3 5])=0;
%% Vecteurs des paramètres connus
%---Paramètres de départ
ini0=zeros(2,nIni)+0.1;
imp0=zeros(2,nImp)+0.1;
inh0=zeros(2,nInh)+0.1;
act0=ones(2,nAct);
par0={ini0 imp0 inh0 act0};
%---Paramètres connus (négatif : paramètres inconnus)
iniKn=zeros(2,nIni);
impKn=zeros(2,nImp);
inhKn=zeros(2,nInh);
actKn=[1;1];
parKn={iniKn impKn inhKn actKn};
%---Paramètres égaux (même numéro=égalité)
%Ini
parApp=0;
parSur=nIni;
iniEq=[parApp+1:parApp+nIni;parSur+1:parSur+nIni];
%Imp
parApp=parSur+nIni;
parSur=parApp+nImp;
impEq=[parApp+1:parApp+nImp;parSur+1:parSur+nImp];
%Inh
parApp=parSur+nImp;
parSur=parApp+nInh;
inhEq=[parApp+1:parApp+nInh;parSur+1:parSur+nInh];
%Act
parApp=parSur+nInh;
parSur=parApp+nAct;
actEq=[parApp+1:parApp+nAct;parSur+1:parSur+nAct];
parEq={iniEq impEq inhEq actEq};
%save('example_CellParam.mat','par0','parKn','parEq')
%% Simulations
%I=I_R;
ini=ini_R;
imp=imp_R;
inh=inh_R;
act=act_R;
%--Paramètres du SBM
XSBM=XSBM_R;
%---Simulation de réseau par SBM
TL=block;
Gr=GSimulSBM(XSBM,TL',k0);
%Valeurs d'initialisation (espèces basales=1, non basales=2)
I=ones(1,n);
I(unique(Gr(Gr(:,3)==1,1)))=2;
%---Etat initial
for i=1:n
    for s=1:Q
        t0(i,s)=rand()<ini(1,I(i));
    end
end
%---Simulation de données
D=dSimul(Gr,A,t0,I,ini,imp,inh,act);
%Calcul de log vraisemblance
lv=lvTot(D,Gr,A,I,ini,imp,inh,act);
%save('example_DataSim.mat','D','Gr','I')
%% Sélection d'arêtes connues : aucune
%Gsel = [i j type(inh ou imp) force present/absent]
Gsel=zeros(0,5);
%% Estimation des paramètres
[iniE,impE,inhE,actE]=parEstKn(D,Gr,A,I,par0,parKn,parEq);
lvEst=lvTot(D,Gr,A,I,ini,imp,inh,act);

%Estimation des paramètres SBM
XSBM0 = [0.01 0.01 0.01 0.01 0.01];
%Changement en matrices d'adjacence
AdjImp=ltoA(Gr(Gr(:,3)==1,[1 2 4]),n,1);
AdjInh=ltoA(Gr(Gr(:,3)==2,[1 2 4]),n,1);
%estimation de alpha
fsbm = @(alpha)-lvSBMalpha(AdjImp, TL, alpha);
alpha = fmincon(fsbm,XSBM0(1), [],[],[],[],0,1,[],optimset('Display','off')); 
%Estimation de Beta1
beta1 = estimbeta1(AdjInh, TL);
if beta1 == 0 beta1 = 0.1; end
%Estimation de Beta2
beta2 = estimbeta2(AdjImp,AdjInh, TL);
if beta2 == 0 beta2 = 0.1; end
%Estimation de alpha et beta fourre-tout
[alphaft, betaft]=estimft(AdjImp, AdjInh, TL);
if alphaft == 0 alphaft = 0.1; end
if betaft == 0 betaft = 0.1; end
%Calcul de la vraisemblance
XSBM_APP = [alpha, beta1, beta2, alphaft, betaft];
lvEst = lvEst +  lvSBM1(AdjImp, AdjInh, TL, XSBM_APP);
%% Restauration
%Ajout de Cplex
addpath(genpath('/opt/ibm/ILOG/CPLEX_Studio1261/cplex/matlab'));
%Lancement de la parallelisation (si pas encore fait)
if isempty(gcp('nocreate'))
    parpool(par);
end
%Sélection d'arcs connus. 
%Pour un arc de i vers j : 
%Gsel = [i j type(inh ou imp) force present/absent]
Gsel=zeros(0,5);
%ILP
lvC=0;
GappCell=cell(n,1);
IappCell=cell(n,1);
%Contraintes de base pour décrire les 
[Ale, Ble, indVar, Ri]=varILP_core(D, k0, nIni, nImp, nInh);
parfor i=1:n
    tic
    GselI=Gsel(Gsel(:,2)==i,:);
    if size(GselI,1)
        Gpres=GselI(GselI(:,5)>0,:);
        Gdeg=Gpres(Gpres(:,5),:);
        k=k0+size(Gdeg,1);
    else
        k=k0
    end
    %Variables de l'ilp (structure de base)
    Ai=Ale;Bi=Ble;
    %Contraintes liées aux arcs connus
    %[A1, B1]=varILP_Gsel(Ai, indVar,GselI);
    %Ai=[Ai; A1];Bi=[Bi B1];
    %Variables de l'ilp (contraintes structurelles)
    [A1, B1]=varILP_structC(Ai, indVar,i, nImp, nInh);
    Ai=[Ai; A1];Bi=[Bi B1];
    %Contraintes liées à la structure d'un réseau écologique
    [A1, B1]=varILP_foodWeb(Ai, indVar);
    Ai=[Ai; A1];Bi=[Bi B1];
    nbV=size(Ai,2);
    %Fonction objectif ILP
    %Ajout des bornes inf et sup comme contraintes
    Ab=speye(nbV);Bb=ones(nbV,1);
    Ab=[Ab;-speye(nbV)];Bb=[Bb;zeros(nbV,1)];
    %Ab=sparse(Ab);
    Ao=[Ai;Ab];Bo=[Bi';Bb];
    Co=lvILP(D, A, i, ini, imp, inh, act, nbV, Ri);
    Co=lvILP_SBM1(D, i, Co, TL', XSBM, indVar);
    [GappI, IappI]=restoration_step_ILP(Ao, Bo, -Co, indVar, i);
    GappCell{i}=GappI;
    IappCell{i}=IappI;
    toc
end
Gapp=cell2mat(GappCell);
Iapp=cell2mat(IappCell);