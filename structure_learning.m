%%%%%%%%%%%%%%%%%%%%%%%
%%%Learning Algorithm%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%Entries%%%
%D:Dataset (n*T*s)
%A:Covariates (n*T*s)
%k:maximum degree
%nodes:Liste de noeuds connus
%ini,imp,inh,act valeur des paramètres
function[G]=structure_learning(D,A,k, ini, imp, inh, act)
n=size(D,1);%nb de noeuds
T=size(D,2);%nb de pas de temps
Q=size(D,3);%nb de jeux de données
%---Param??tres du DBN (utilises pour simuler les donnees)
ini=[0.8 0.8;1 0];
imp=[0 0;0.8 0.8];
inh=[0;0.8];
act=[0.8 1;0.8 1];
%---Etat initial
for i=1:n
    t0(i)=rand()<ini(1,I(i));
end
%---Réseau (sous forme de liste : n1 n2 type(1=imp 2=inh) etiquette
G=[ 1 2 1 1;
    1 3 1 2;
    1 4 1 1;
    2 3 1 1;
    2 1 2 1;
    3 2 2 1;
    4 1 2 1];
%---Simulation de données
D=dSimul(G,A,t0,I,ini,imp,inh,act);
%Calcul de log vraisemblance
lv=lvTot(D,G,A,I,ini,imp,inh,act);
%Estimation des paramètres
%---Paramètres de départ
ini0=[0.01 0.01;0.01 0.01];
imp0=[0.01 0.01;0.01 0.01];
inh0=[0.01;0.01];
act0=[0.01 0.01;0.01 0.01];
[iniE,impE,inhE,actE]=parEst(D,G,A,I,ini0,imp0,inh0,act0);
% Données non simulées
D = [0     0     1     1     1     1     1;
     0     0     1     0     1     1     1;
     0     1     0     1     1     0     1;
     0     0     0     0     1     0     1];
%Calcul de log vraisemblance
lv=lvTot(D,G,A,I,ini,imp,inh,act);
 %Variables nécéssaires pour la construction de l'ILP
Gsel=[G(1,:) 1];
nImp=2;
nInh=1;
nIni=2;
%Variables ILP
lvC=0;
for i=1:n
    GselI=Gsel(find(Gsel(:,2)==i),:);
    Gpres=GselI(find(GselI(:,5)>0),:);
    Gdeg=Gpres(find(Gpres(:,5)),:);
    k=k0+size(Gdeg,1);
    %Variables de l'ilp (structure de base)
    [Ale, Ble, indVar, Ri]=varILP_core(D, k, nIni, nImp, nInh);
    %Variables de l'ilp (contraintes structurelles)
    [A1, B1]=varILP_structC(Ale, indVar,i, nImp, nInh);
    Ale=[Ale; A1];Ble=[Ble B1];
    %Variables de l'ilp (interaction connues
    [A1, B1]=varILP_Gsel(Ale, indVar, GselI);
    Ale=[Ale; A1];Ble=[Ble B1];
    %Variables de l'ilp (interaction connues
    [A1, B1]=varILP_foodWeb(Ale,indVar);
    Ale=[Ale; A1];Ble=[Ble B1];
    nbV=size(Ale,2);
    %Fonction objectif ILP
    %Ajout des bornes inf et sup comme contraintes
    Ab=speye(nbV);Bb=ones(nbV,1);
    Ab=[Ab;-speye(nbV)];Bb=[Bb;zeros(nbV,1)];
    %Ab=sparse(Ab);
    Ao=[Ale;Ab];Bo=[Ble';Bb];
    C=lvILP(D, A, i, ini, imp, inh, act, nbV, Ri);
    Gi=G(find(G(:,2)==i),:);
    Ax=XvarILP(D, Gi, GselI,i, I(i), k, nbV, nIni, nImp, nInh);
    lvC=lvC+(C*Ax');
    %Gapp=restoration_step(D,A,k,i, ini, imp, inh,act);
end
lv
lvC