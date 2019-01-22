%varILP_Gsel(A0, indVar, GselI) :
%Génère les contrainte de l'ILP ajoutant un ensemble d'intéractions connues
%% INPUT
% A0 : Matrice A de l'ILP généré par la fonction varILP_core
% indVar : Indices des variables de l'ILP généré par la fonction varILP_core
% GselI : Vecteur renseignant les itnéractions connues de la forme : [noeud1, noeud2, type(1=impulseur 2=inhibiteur), étiquette, presence(1=presence de l'arête, 0=absence de l'arête)]
%% OUTPUT
% A1 : Matrice A
% B1 : Vecteur B
function[A1, B1]=varILP_Gsel(A0, indVar, GselI)
    A1=sparse(0,size(A0,2));
    B1=[];
    Gimp=indVar{1,1};
    Ginh=indVar{1,2};
    rowC=0;%Indice de la derni??re contrainte (la ligne ?? ajouter dans A et B), colC est la colonne
        %Contrainte du graphe connu : les liens connus sont mis en contrainte
        %Gif
        for j=1:size(GselI,1)
            if GselI(j,5)
                if GselI(j,3)==1
                    rowC=rowC+1;
                    A1(rowC,Gimp(GselI(j,4),GselI(j,1)))=-1;
                    B1(rowC)=-1;
                end
                if GselI(j,3)==2
                    rowC=rowC+1;
                    A1(rowC,Ginh(GselI(j,4), GselI(j,1)))=-1;
                    B1(rowC)=-1;
                end
            else
                if GselI(j,3)==1
                    rowC=rowC+1;
                    A1(rowC,Gimp(GselI(j,4), GselI(j,1)))=1;
                    B1(rowC)=0;
                end
                if GselI(j,3)==2
                    rowC=rowC+1;
                    A1(rowC,Ginh(GselI(j,4), GselI(j,1)))=1;
                    B1(rowC)=0;
                end
            end
        end