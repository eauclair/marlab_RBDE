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
function[Ax]=XvarILP(D, Gi, GselI, i, Ui, k, nbV, nIni, nImp, nInh)
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
    for imp=1:nImp
        temp=[1:n]+last;%indices des inhibiteurs (i,j)
        last=last+n;
        Gimp=[Gimp;temp];
    end
    Ginh=[];
    for imp=1:nInh
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
        for imp=1:nImp
            Mimp(:,:,d,imp)=reshape(1:T*K,[T,K])+last;%indices de Mf(i,t,d)
            last=last+T*K;
            Nimp(:,:,d,imp)=reshape(1:T*K,[T,K])+last;%indices de Nf(i,t,d)
            last=last+T*K;       
            Limp(:,:,d,imp)=reshape(1:T*K,[T,K])+last;%indices de Lf(i,t,d)
            last=last+T*K;
        end
        for inh=1:nInh
            Minh(:,:,d,inh)=reshape(1:T*K,[T,K])+last;%indices de M-(i,t,d)
            last=last+T*K;
            Ninh(:,:,d,inh)=reshape(1:T*K,[T,K])+last;%indices de N-(i,t,d)
            last=last+T*K;
            Linh(:,:,d,inh)=reshape(1:T*K,[T,K])+last;%indices de L-(i,t,d)
            last=last+T*K;
        end
    end
    %Nombre de variables
    %Pr??allocation : matrice A
    Ax=sparse(zeros(1,nbV));
    %Pr??allocation : matrice B
    Ble=zeros(1,1);
    %Pr??allocation : Matrice pour s'y retrouver dans les R
    %Ri=zeros((factoriaVestl(k+3)/(6*factorial(k))+2*(factorial(k+2)/(2*factorial(k)))+k+1)*T,7);
    Ri=[];
    rowC=0;%Indice de la derni??re contrainte (la ligne ?? ajouter dans A et B), colC est la colonne
    ir=0;%Indice de la matrice Ri
    for j=1:size(Gi,1)
        if Gi(j,3)==1
            Ax(Gimp(Gi(j,4),Gi(j,1)))=1;
        else
            Ax(Ginh(Gi(j,4),Gi(j,1)))=1;
        end
    end
    Ax(Init(Ui))=1;
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
                Gjimp=Gi(find(Gi(:,3)==1),:);
                Gjinh=Gi(find(Gi(:,3)==2),:);
                for imp=1:nImp
                    %Mimp*(d+1)<=sum(gfji*xjt)+1
                    Gj=Gjimp(find(Gjimp(:,4)==imp),:);
                    Dj=D(Gj(:,1),t,s);
                    if sum(Dj)>=dImp(imp)
                        Ax(Mimp(t,dImp(imp)+1,s,imp))=1;
                    end
                    if sum(Dj)<=dImp(imp)
                        Ax(Nimp(t,dImp(imp)+1,s,imp))=1;
                    end
                    if sum(Dj)==dImp(imp)
                        Ax(Limp(t,dImp(imp)+1,s,imp))=1;
                    end
                end
                for inh=1:nInh
                    %Minh*(d+1)<=sum(gfji*xjt)+1
                    Gj=Gjinh(find(Gjinh(:,4)==inh),:);
                    Dj=D(Gj(:,1),t,s);
                    if sum(Dj)>=dInh(inh)
                        Ax(Minh(t,dInh(inh)+1,s,inh))=1;
                    end
                    if sum(Dj)<=dInh(inh)
                        Ax(Ninh(t,dInh(inh)+1,s,inh))=1;
                    end
                    if sum(Dj)==dInh(inh)
                        Ax(Linh(t,dInh(inh)+1,s,inh))=1;
                    end
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
                    if Ui==ini
                        %Nombre d'étiquettes ayant le bon nombre de parents
                        %présents
                        nbLab=0;
                        %R<=Limp(dimp)
                        for imp=1:nImp
                            if Ax(Limp(t,dImp(imp)+1,s,imp))
                                nbLab=nbLab+1;
                            end
                        end
                        %Limp(t,dImp(imp)+1,s,imp)
                        %R<=Linh(dinh)
                        for inh=1:nInh
                            if Ax(Linh(t,dInh(inh)+1,s,inh))
                                %Linh(t,dInh(inh)+1,s,inh)
                                %Ax(Ri(ir,1))=1;
                                nbLab=nbLab+1;
                            end
                        end
                        if nbLab==nImp+nInh
                            Ax(Ri(ir,1))=1;
                        end
                    end
                end
                %Augmentation d'une valeur d'une certaine étiquette : 
                %On augmente la valeur d d'une seule étiquette, et on sort
                %de la boucle si on ne peut plus augmenter l'étiquette sans
                %dépasser k.
                %Augmentation d'une valeur d'une certaine étiquette : 
                %On augmente la valeur d d'une seule étiquette, et on sort
                %de la boucle si on ne peut plus augmenter l'étiquette sans
                %dépasser k.
%                 modif=0;                    
%                 label=1;
%                 while ~modif
%                     if sum(dImp)+sum(dInh)<K-1
%                         dImp(label)=dImp(label)+1;
%                         modif=1;
%                     else
%                         label=label+1;
%                         dImp(label-1)=0;
%                         if label>nImp
%                             label=1;
%                             if sum(dImp)+sum(dInh)<K-1
%                                 dInh(label)=dInh(label)+1;
%                                 modif=1;
%                             else
%                                 label=label+1;
%                                 if label>nInh
%                                     modif=1;
%                                     go=0;
%                                 else
%                                     dInh(label-1)=0;
%                                 end
%                             end
%                         else
%                             dImp(label-1)=0;
%                         end
%                     end
%                 end
            %Fin de la boucle pour augmenter la valeur d'une étiquette
                nPar=nImp+nInh;
                dPar=[dImp dInh];
                modif=0; 
                label=1;
                while ~modif
                    if sum(dPar)<K-1
                        dPar(label)=dPar(label)+1;
                        modif=1;
                    else
                        label=label+1;
                        if label>nPar
                            modif=1;
                            go=0;
                        else
                          dPar(label-1)=0;
                        end
                    end
                end
                if go
                    dImp=dPar(1:nImp);
                    dInh=dPar(nImp+1:nImp+nInh);
                end
            %Fin de la boucle pour augmenter la valeur d'une étiquette
            end
        end
    end 
%A chaque esp??ce et chaque pas de temps doit ??tre associ?? une variable R
%rowC=rowC+1;
%Ale(rowC, last+1:last+ir)=-1;
%Ble(rowC)=-n*T;