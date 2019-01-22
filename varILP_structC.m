%varILP_structC(A0, indVar, GselI) :
%Génère les contrainte de l'ILP ajoutant un ensemble de contraintes sur le
%graphe très classique : une seule intéraction entre 2 espèces, et aucune
%interaction entre une espèce et elle-même
%% INPUT
% A0 : Matrice A de l'ILP généré par la fonction varILP_core
% indVar : Indices des variables de l'ILP généré par la fonction varILP_core
% i : indice de l'espèce dont on cherche à apprendre les parents
% nImp : Nombre d'impulseurs dans le modèle
% nInh : Nombre d'inhibiteurs dans le modèle
%% OUTPUT
% A1 : Matrice A
% B1 : Vecteur B
function[A1, B1]=varILP_structC(A0,indVar, i, nImp, nInh)
    A1=sparse(0,size(A0,2));
    B1=[];
    Gimp=indVar{1,1};
    Ginh=indVar{1,2};
    n=size(Gimp,2);
    rowC=0;%Indice de la derni??re contrainte (la ligne ?? ajouter dans A et B), colC est la colonne
        %Un seul lien entre deux esp??ces au maximum
        %sum(Gjiimp)+sum(Gjiinh)<=1
        for j=1:n
            rowC=rowC+1;
            A1(rowC,Gimp(:,j))=1;
            A1(rowC,Ginh(:,j))=1;
            B1(rowC)=1;
        end    
        %Pas de lien r??ccursif (self loop)
        %Giip+Giif+Gii-<=0
        rowC=rowC+1;
        for idImp=1:nImp
            A1(rowC,Gimp(idImp,i))=1;
        end
        for idInh=1:nInh
            A1(rowC,Ginh(idInh,i))=1;
        end
        B1(rowC)=0;