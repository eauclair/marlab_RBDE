%% Apprentissage R-E avec connaissances brutes
%Demande : k0 : degré max ;   D : Données présence/absence ; 
%I : Initialisation ; par0, parEq, parKn : paramètres ;
% SBMBlock : bloc du SBM ; Gsel : connaissances "dures" des arcs connus
% (peut être vide)
n=size(D,1);
T=size(D,2);
Q=size(D,3);
%Graphe vide
G=zeros(0,4);
[iniE,impE,inhE,actE]=parEstKn(D,G,A,I,par0,parKn,parEq);
lv0=lvTot(D,G,A,I,iniE,impE,inhE,actE);
%% Apprentissage
%Ajout de Cplex
addpath(genpath('/opt/ibm/ILOG/CPLEX_Studio1261/cplex/matlab'));
%Lancement de la parallelisation (si pas encore fait)
if isempty(gcp('nocreate'))
    parpool(par);
end
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
    %Estimation des paramètres
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