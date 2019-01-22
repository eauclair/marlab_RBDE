%% Apprentissage R-E avec connaissances brutes
function[graphFinTot, IniFinTot, timeERFinTot, lvFinTot, paramFinTot]=ERlearning(D,A,I,k0,par0,parKn,parEq,Gsel)
n=size(D,1);
nIni=size(par0{1},2);
nImp=size(par0{2},2);
nInh=size(par0{3},2);
%Graphe vide
G=Gsel(Gsel(:,5)==1,1:4);
%Première étape E
[iniE,impE,inhE,actE]=parEstKn(D,G,A,I,par0,parKn,parEq);
lv0=lvTot(D,G,A,I,iniE,impE,inhE,actE);
go=1;
timeER=tic;
while go
    %-------------Etape R-----------------
    lvC=0;
    GappCell=cell(n,1);
    IappCell=cell(n,1);
    %Contraintes de base pour décrire les variables ILP
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
        nbV=size(Ai,2);
        %Fonction objectif ILP
        %Ajout des bornes inf et sup comme contraintes
        Ab=speye(nbV);Bb=ones(nbV,1);
        Ab=[Ab;-speye(nbV)];Bb=[Bb;zeros(nbV,1)];
        Ao=[Ai;Ab];Bo=[Bi';Bb];
        Co=lvILP(D, A, i, iniE, impE, inhE, actE, nbV, Ri);
        [GappI, IappI]=restoration_step_ILP(Ao, Bo, -Co, indVar, i);
        GappCell{i}=GappI;
        IappCell{i}=IappI;
    end
    Gapp=cell2mat(GappCell);
    Iapp=cell2mat(IappCell);
    %-------------Etape E-----------------
    [iniE,impE,inhE,actE]=parEstKn(D,Gapp,A,Iapp,par0,parKn,parEq);
    %Calcul de la vraisemblance
    lv1=lvTot(D,Gapp,A,Iapp,iniE,impE,inhE,actE);
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