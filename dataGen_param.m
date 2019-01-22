%% Fonction de génération de données : paramètres
function[ini, imp, inh, act, par0, parKn, parEq]=dataGen_param(nIni, nImp, nInh, nAct, parExist)
val0=0.001;
%% Initialisation des variables
%Valeurs de paramètres
ini=zeros(2,nIni);
imp=zeros(2,nImp);
inh=zeros(2,nInh);
act=ones(2,nAct);
%Valeur initiales
ini0=zeros(2,nIni);
imp0=zeros(2,nImp);
inh0=zeros(2,nInh);
act0=ones(2,nAct);
%Paramètres connus
iniKn=ones(2,nIni);
impKn=ones(2,nImp);
inhKn=ones(2,nInh);
actKn=ones(2,nAct);
%% Valeurs de paramètres aléatoires
if parExist(1,1)
    ini(1,:)=rand(1,nIni);
    ini0(1,:)=zeros(1,nIni)+val0;
    iniKn(1,:)=zeros(1,nIni);
end
if parExist(2,1)
    ini(2,:)=rand(1,nIni);
    ini0(2,:)=zeros(1,nIni)+val0;
    iniKn(2,:)=zeros(1,nIni);
end
if parExist(1,2)
    imp(1,:)=rand(1,nImp);
    imp0(1,:)=zeros(1,nImp)+val0;
    impKn(1,:)=zeros(1,nImp);
end
if parExist(2,2)
    imp(2,:)=rand(1,nImp);
    imp0(2,:)=zeros(1,nImp)+val0;
    impKn(2,:)=zeros(1,nImp);
end
if parExist(1,3)
    inh(1,:)=rand(1,nInh);
    inh0(1,:)=zeros(1,nInh)+val0;
    inhKn(1,:)=zeros(1,nInh);
end
if parExist(2,3)
    inh(2,:)=rand(1,nInh);
    inh0(2,:)=zeros(1,nInh)+val0;
    inhKn(2,:)=zeros(1,nInh);
end
if parExist(1,4)
    act(1,:)=rand(1,nAct);
    act0(1,:)=zeros(1,nAct)+val0;
    actKn(1,:)=zeros(1,nAct);
end
if parExist(2,4)
    act(2,:)=rand(1,nAct);
    act0(2,:)=zeros(1,nAct)+val0;
    actKn(2,:)=zeros(1,nAct);
end
par0={ini0 imp0 inh0 act0};
parKn={iniKn impKn inhKn actKn};
%% Paramètres égaux (même numéro=égalité)
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