%varOptGSBMEsp(D, G, Gsel, A, i, k, TL) :
%genere les variables et les contrainte de l'ILP d'une
%esp?ce (A.x <= B), dans le cas ou on a de la connaissance
% sur la presence ou l'absence de certaines aretes.

%% INPUT
% D : donnees,
% Gr : le vrai voisinage de l'esp�ce i, taille nx3.
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
function[Ale, Ble, Ri, Gimp, Ginh,Init]=varILP(D, GselI, i, k, nIni, nImp, nInh)
    n=size(D,1);
    T=size(D,2);
    Q=size(D,3);
%     Gpres=GselI(find(GselI(:,5)>0),:);
%     Gdeg=Gpres(find(Gpres(:,5)),:);
%     K=k+size(Gdeg,1)+1;
    K=k+1;
    last=0;
    %Gp=1:n+last;%indices de Gp(i,j)
    %Ce sont ici des lignes, il faut des colonnes dans les matrices de
    %contingence
    %last=last+n;
    Gimp=[];
    for idImp=1:nImp
        temp=[1:n]+last;%indices des inhibiteurs (i,j)
        last=last+n;
        Gimp=[Gimp;temp];
    end
    Ginh=[];
    for idInh=1:nInh
        temp=[1:n]+last;%indices de Gimp(imp,:)(i,j)
        last=last+n;
        Ginh=[Ginh;temp];
    end 
    %Gimpinh = [1:n] + last; % indice des Ginhimp, egal � 1 ssi il y a à la fois un inhibiteur et un impulseur de j vers i
    %last = last +n;
    %Liste des initialisations possibles
    Init=[1:nIni]+last;
    last=last+nIni;
    %Matrices des nombres de parents (borne inf/borne sup/exact)
    Mimp=[];
    Minh=[];
    Nimp=[];
    Ninh=[];
    Limp=[];
    Linh=[];
    for d=1:Q
        for idImp=1:nImp
            Mimp(:,:,d,idImp)=reshape(1:T*K,[T,K])+last;%indices de Mf(t,d,s,lImp)
            last=last+T*K;
            Nimp(:,:,d,idImp)=reshape(1:T*K,[T,K])+last;%indices de Nf(t,d,s,lImp)
            last=last+T*K;       
            Limp(:,:,d,idImp)=reshape(1:T*K,[T,K])+last;%indices de Lf(t,d,s,lImp)
            last=last+T*K;
        end
        for idInh=1:nInh
            Minh(:,:,d,idInh)=reshape(1:T*K,[T,K])+last;%indices de M-(t,d,s,lInh)
            last=last+T*K;
            Ninh(:,:,d,idInh)=reshape(1:T*K,[T,K])+last;%indices de N-(t,d,s,lInh)
            last=last+T*K;
            Linh(:,:,d,idInh)=reshape(1:T*K,[T,K])+last;%indices de L-(t,d,s,lInh)
            last=last+T*K;
        end
    end
    %Nombre de variables
    nbV=last;%+(factorial(k+3)/(6*factorial(k))+2*(factorial(k+2)/(2*factorial(k)))+k+1)*n*T;
    %Pr??allocation : matrice A
    Ale=sparse(1, nbV);
    %Pr??allocation : matrice B
    Ble=zeros(1,1);
    %Pr??allocation : Matrice pour s'y retrouver dans les R
    %Ri=zeros((factoriaVestl(k+3)/(6*factorial(k))+2*(factorial(k+2)/(2*factorial(k)))+k+1)*T,7);
    Ri=[];
    rowC=0;%Indice de la derni??re contrainte (la ligne ?? ajouter dans A et B), colC est la colonne
    ir=0;%Indice de la matrice Ri
        %Quantit??s d??pendantes uniquement de i
        %Le maximum de parents est k
        %sum(Gjip+Gjif+Gji-)<=k
        rowC=rowC+1;
        Ale(rowC,[Gimp(:)',Ginh(:)'])=1;
        Ble(rowC)=k;
        %Un seul lien entre deux esp??ces au maximum
        %sum(Gjiimp)+sum(Gjiinh)<=1
        for j=1:n
            rowC=rowC+1;
            Ale(rowC,Gimp(:,j))=1;
            Ale(rowC,Ginh(:,j))=1;
            Ble(rowC)=1;
        end
        %Contrainte du graphe connu : les liens connus sont mis en contrainte
        %Gif
        for j=1:size(GselI,1)
            if GselI(j,5)
                if GselI(j,3)==1
                    rowC=rowC+1;
                    Ale(rowC,Gimp(GselI(j,1),GselI(j,4)))=-1;
                    Ble(rowC)=-1;
                end
                if GselI(j,3)==2
                    rowC=rowC+1;
                    Ale(rowC,Ginh(GselI(j,1),GselI(j,4)))=-1;
                    Ble(rowC)=-1;
                end
            else
                if GselI(j,3)==1
                    rowC=rowC+1;
                    Ale(rowC,Gimp(GselI(j,1),GselI(j,4)))=1;
                    Ble(rowC)=0;
                end
                if GselI(j,3)==2
                    rowC=rowC+1;
                    Ale(rowC,Ginh(GselI(j,1),GselI(j,4)))=1;
                    Ble(rowC)=0;
                end
            end
        end
        %Pas de lien r??ccursif (self loop)
        %Giip+Giif+Gii-<=0
        rowC=rowC+1;
        for idImp=1:nImp
            Ale(rowC,Gimp(idImp,i))=1;
        end
        for idInh=1:nInh
            Ale(rowC,Ginh(idInh,i))=1;
        end
        Ble(rowC)=0;
        %Une seule initialisation
        %Giip+Giif+Gii-<=0
        rowC=rowC+1;
        for ini=1:nIni
            Ale(rowC,Init(ini))=1;
        end
        Ble(rowC)=1;
    for s=1:Q
        %-------------------------------
        for t=1:T-1
            %Quantit??s d??pendantes de t et de i
            %Quantit??s d??pendantes du nombre de parents de chaque type
            dImp=zeros(1,nImp);
            dInh=zeros(1,nInh);
            %test=[];%Juste pour l'affichage
            go=1;
            while go
                %test=[test;dImp dInh];
                for idImp=1:nImp
                    %Mimp*(d+1)<=sum(gfji*xjt)+1
                    rowC=rowC+1;
                    Ale(rowC, Mimp(t,dImp(idImp)+1,s,idImp))=dImp(idImp)+1;
                    Ale(rowC,Gimp(idImp,:))=-D(:,t,s)';
                    Ble(rowC)=1;
                    %Mimp*(K-d)>sum(gfji*xjt)-d
                    rowC=rowC+1;
                    Ale(rowC, Mimp(t,dImp(idImp)+1,s,idImp))=-K+dImp(idImp);
                    Ale(rowC,Gimp(idImp,:))=D(:,t,s)';
                    Ble(rowC)=dImp(idImp)-1;
                    %--Nimp--
                         %Nimp*(K-d)<=K-sum(gfji*xjt)
                    rowC=rowC+1;
                    Ale(rowC, Nimp(t,dImp(idImp)+1,s,idImp))=K-dImp(idImp);
                    Ale(rowC,Gimp(idImp,:))=D(:,t,s)';
                    Ble(rowC)=K;
                        %Nimp*(d+1)>d-sum(gfji*xjt)
                    rowC=rowC+1;
                    Ale(rowC, Nimp(t,dImp(idImp)+1,s,idImp))=-dImp(idImp)-1;
                    Ale(rowC,Gimp(idImp,:))=-D(:,t,s)';
                    Ble(rowC)=-dImp(idImp)-1;
                    %--Limp--
                        %Limp<=Mimp
                    rowC=rowC+1;
                    Ale(rowC, Limp(t,dImp(idImp)+1,s,idImp))=1;
                    Ale(rowC,Mimp(t,dImp(idImp)+1,s,idImp))=-1;
                    Ble(rowC)=0;
                        %Limp<=Nimp
                    rowC=rowC+1;
                    Ale(rowC, Limp(t,dImp(idImp)+1,s,idImp))=1;
                    Ale(rowC,Nimp(t,dImp(idImp)+1,s,idImp))=-1;
                    Ble(rowC)=0;
                        %Lf>=Mf+Nf-1
                    rowC=rowC+1;
                    Ale(rowC, Limp(t,dImp(idImp)+1,s,idImp))=-1;
                    Ale(rowC,Mimp(t,dImp(idImp)+1,s,idImp))=1;
                    Ale(rowC,Nimp(t,dImp(idImp)+1,s,idImp))=1;
                    Ble(rowC)=1;
                end
                for idInh=1:nInh
                   %--Mn--
                        %Mn*(d+1)<=sum(gnji*xjt)+1
                    rowC=rowC+1;
                    Ale(rowC, Minh(t,dInh(idInh)+1,s,idInh))=dInh(idInh)+1;
                    Ale(rowC,Ginh(idInh,:))=-D(:,t,s)';
                    Ble(rowC)=1;
                        %Mn*(K-d)>sum(gnji*xjt)-d
                    rowC=rowC+1;
                    Ale(rowC, Minh(t,dInh(idInh)+1,s))=-K+dInh(idInh);
                    Ale(rowC,Ginh(idInh,:))=D(:,t,s)';
                    Ble(rowC)=dInh(idInh)-1;
                    %--Nn--
                         %Nn*(K-d)<=K-sum(gnji*xjt)
                    rowC=rowC+1;
                    Ale(rowC, Ninh(t,dInh(idInh)+1,s))=K-dInh(idInh);
                    Ale(rowC,Ginh(idInh,:))=D(:,t,s)';
                    Ble(rowC)=K;
                        %Nn*(d+1)>d-sum(gnji*xjt)
                    rowC=rowC+1;
                    Ale(rowC, Ninh(t,dInh(idInh)+1,s))=-dInh(idInh)-1;
                    Ale(rowC,Ginh(idInh,:))=-D(:,t,s)';
                    Ble(rowC)=-dInh(idInh)-1;
                    %--Ln--
                        %Ln<=Mn
                    rowC=rowC+1;
                    Ale(rowC, Linh(t,dInh(idInh)+1,s))=1;
                    Ale(rowC,Minh(t,dInh(idInh)+1,s))=-1;
                    Ble(rowC)=0;
                        %Ln<=Nn
                    rowC=rowC+1;
                    Ale(rowC, Linh(t,dInh(idInh)+1,s))=1;
                    Ale(rowC,Ninh(t,dInh(idInh)+1,s))=-1;
                    Ble(rowC)=0;
                        %Ln>=Mn+Nn-1
                    rowC=rowC+1;
                    Ale(rowC, Linh(t,dInh(idInh)+1,s))=-1;
                    Ale(rowC,Minh(t,dInh(idInh)+1,s))=1;
                    Ale(rowC,Ninh(t,dInh(idInh)+1,s))=1;
                    Ble(rowC)=1;
                end
                for ini=1:nIni
                    %--R--
                    ir=ir+1;
                    Ri(ir, 1)=last+ir;%Indice
                    Ri(ir, 2)=t;%Temps
                    Ri(ir,3)=s;%Jeu de données
                    Ri(ir,4)=ini;%Type d'initialisation
                    ind(1)=4;
                    Ri(ir, 1+ind(1):nImp+ind(1))=dImp;%Nb de facilitateurs vivants
                    ind(2)=ind(1)+nImp;
                    Ri(ir, 1+ind(2):nInh+ind(2))=dInh;%Nb de facilitateurs vivants
                    %R<=ini
                    rowC=rowC+1;
                    Ale(rowC,Ri(ir,1))=1;
                    Ale(rowC,Init(ini))=-1;
                    Ble(rowC)=0;
                    %R<=Limp(dimp)
                    for idImp=1:nImp
                        rowC=rowC+1;
                        Ale(rowC,Ri(ir,1))=1;
                        Ale(rowC,Limp(t,dImp(idImp)+1,s,idImp))=-1;
                        Ble(rowC)=0;
                    end
                    %R<=Linh(dinh)
                    for idInh=1:nInh
                        rowC=rowC+1;
                        Ale(rowC,Ri(ir,1))=1;
                        Ale(rowC,Linh(t,dInh(idInh)+1,s,idInh))=-1;
                        Ble(rowC)=0;
                    end
                    %R>=ini+Limp(dimp)+Linh(dinh)-(nImp+nInh+1)
                    rowC=rowC+1;
                    Ale(rowC,Ri(ir,1))=-1;
                    Ale(rowC,Init(ini))=1;
                    for idImp=1:nImp
                        Ale(rowC,Limp(t,dImp(idImp)+1,s,idImp))=1;
                    end
                    for idInh=1:nInh
                        Ale(rowC,Linh(t,dInh(idInh)+1,s,idInh))=1;
                    end
                    %Nombre de variables -1, soit nImp+nInh+1-1
                    Ble(rowC)=nImp+nInh;
                end
                %Augmentation d'une valeur d'une certaine étiquette : 
                %On augmente la valeur d d'une seule étiquette, et on sort
                %de la boucle si on ne peut plus augmenter l'étiquette sans
                %dépasser k.
                modif=0;                    
                label=1;
                while ~modif
                    if sum(dImp)+sum(dInh)<K-1
                        dImp(label)=dImp(label)+1;
                        modif=1;
                    else
                        label=label+1;
                        dImp(label-1)=0;
                        if label>nImp
                            label=1;
                            if sum(dImp)+sum(dInh)<K-1
                                dInh(label)=dInh(label)+1;
                                modif=1;
                            else
                                label=label+1;
                                if label>nInh
                                    modif=1;
                                    go=0;
                                else
                                    dInh(label-1)=0;
                                end
                            end
                        else
                            dImp(label-1)=0;
                        end
                    end
                end
            %Fin de la boucle pour augmenter la valeur d'une étiquette
            end
        end
    end 
%A chaque esp??ce et chaque pas de temps doit ??tre associ?? une variable R
%rowC=rowC+1;
%Ale(rowC, last+1:last+ir)=-1;
%Ble(rowC)=-n*T;