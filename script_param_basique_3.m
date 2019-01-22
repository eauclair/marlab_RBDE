%% Script de génération de paramètres et d'informations complémantaire
%1 : epsilon pour l'apparition ; rho et tau pour la survie
%2 : rho et tau pour l'apparition ; epsilon pour la survie
%3 : epsilon, rho et tau pour l'apparition et la survie
%% Covariables, blocs et comportements inhérents aléatoires
%Covariable
A_R=unidrnd(nAct,n,T,Q);
%Types d'initialisation
I_R=unidrnd(nIni,n,1);
%Blocs aléatoires
Block_R=unidrnd(nBlock, n, 1);
% Paramètres SBM
XSBM_F=[0.8 0.3 0.2 0.4 0.4];
%% Paramètres générés aléatoirement
%---Param??tres du RBDE
ini_R=rand(2,nIni);
imp_R=rand(2,nImp);
inh_R=rand(2,nInh);
act_R=ones(2,nAct);
%% Vecteurs des paramètres connus
%---Paramètres de départ
ini0=zeros(2,nIni)+0.1;
imp0=zeros(2,nImp)+0.1;
inh0=zeros(2,nInh)+0.1;
act0=ones(2,nAct);
%act0=zeros(2,nAct)+0.1;
par0={ini0 imp0 inh0 act0};
%---Paramètres connus (négatif : paramètres inconnus)
iniKn=zeros(2,nIni);
impKn=zeros(2,nImp);
inhKn=zeros(2,nInh);
actKn=ones(2,nAct);
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
%Paramètres
ini=ini_R;
imp=imp_R;
inh=inh_R;
act=act_R;
%Paramètres du SBM
XSBM=XSBM_F;
%Blocs
SBMblock=Block_R;
%Covariables
A=A_R;