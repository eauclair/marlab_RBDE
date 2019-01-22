%% Script de simulation de graphe et de données
%Demande : 
%---XSBM : paramètres du SBM
%---SBMblock : blocks du SBM
%---k0 : degré maximal
%---n, T, Q : Nombre de noeuds, de pas de temps et de JDD indépendants
%---A : covariables de chaque noeud
%---I : Initialisation de chaque noeud
%---ini, imp, inh, act : valeurs des paramètres RBDE

%Graphe simulé
Gr=GSimulSBM(XSBM,SBMblock',k0);
%Valeurs d'initialisation (aléatoires)
I=I_R;
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