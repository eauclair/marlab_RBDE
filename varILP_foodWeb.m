%varOptGSBMEsp(D, G, Gsel, A, i, k, TL) :
%genere les variables et les contrainte de l'ILP d'une
%esp?ce (A.x <= B), dans le cas ou on a de la connaissance
% sur la presence ou l'absence de certaines aretes.

%% INPUT
% D : donnees,
% Gr : le vrai voisinage de l'espece i, taille nx3.
% Gsel : vecteur binaire colonne (nx1) , avec 1 en j si on connait la
%        valeur de gji et 0 sinon.
% A : covariable,
% i : indice del'espece
% k : nombre de parents max
% TL : vecteur  colonne des niveaux trophiques de chaque espece

%% OUTPUT
% Ale : Matrice A contenant les quantites variables des contraintes
% Ble : Matrice B contenant la valeur des contraintes
% Ri :  descrition des variables Ri
                    %Ri(ir, 1) Temps
                    %Ri(ir, 2) Nb de proies vivantes
                    %Ri(ir, 3) Nb de facilitateurs vivants
                    %Ri(ir, 4) Nb de n??gatifs vivants
                    %Ri(ir, 5) Comportement trophique
                    %Ri(ir, 6) indice de cette variable;
% Gi : matrice 3 by n avec Gi(:,2) = indices des variables gjif et
% Gi(:,3) indice des variables gjin 
% Gfn : vecteur 1 by n avec les indices des variables gjifn

%%!!! ATTENTION !!!!!!!!!!!!!!!!!
%Cette fonction n'apprendra que des relations de type facilitation et
%influence n??gative (+ et -)

% NP 23/01/17 correction bug formule dans la partie SBM
function[A1, B1]=varILP_foodWeb(A0,indVar)
    A1=sparse(0,size(A0,2));
    B1=[];
    Gimp=indVar{1,1};
    Init=indVar{1,3};
    n=size(Gimp,2);
    rowC=0;%Indice de la derni??re contrainte (la ligne ?? ajouter dans A et B), colC est la colonne 
    %I(2)<=sum(Gimp)
    rowC=rowC+1;
    A1(rowC,Init(2))=1;
    A1(rowC,Gimp)=-1;
    B1(rowC)=0;
    %I(2)>=Gimp(1),Gimp(2)...
    for j=1:n
        rowC=rowC+1;
        A1(rowC,Init(2))=-1;
        A1(rowC,Gimp(:,j))=1;
        B1(rowC)=0;
    end