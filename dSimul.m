%Dsimul(G,X, A)
%Simule des données d'observation à partir d'un graphe G à l'aide des lois
%de probas définies par les paramètres ini, imp, inh, act
%A décrit l'action sur chaque noeud à chaque pas de temps, et I le type d'initialisation
%t0 : vecteur de taille n*s : état initial à t=0
function[D] = dSimul(G,A,t0,I,ini,imp,inh,act)
    %variables
    n=size(A,1);
    T=size(A,2);
    Q=size(A,3);
    D=zeros(n, T+1, Q);
    D(:,1,:)=t0;
    for i=1:n
        Gi=G(find(G(:,2)==i),:);
        for t=1:T
            for s=1:Q
                %Absence au départ
                if D(i,t,s)==0
                    Pini=ini(1,I(i));%Apparition spontanée (On ne donne en entrée que les paramètres d'initialisation relatifs au noeud i)
                    Pimp=1;%Initialisation de la proba (impulseurs)
                    Pinh=1;%Initialisation de la proba (inhibiteurs)
                    Pact=act(1,A(i,t,s));%Proba relative à l'action sur i au temps t dans le jeu de données s
                    for j=1:size(Gi,1)%j est l'identification d'une intéraction sur i
                        if Gi(j,3)==1%Cette intéraction est un impulseur
                            if D(Gi(j,1),t,s)%Le noeud impulseur est présent
                                Pimp=Pimp*(1-imp(1,Gi(j,4)));%Probabilité que l'impulseur j>i échoue
                            end
                        end
                        if Gi(j,3)==2%Cette intéraction est un inhibiteur
                            if D(Gi(j,1),t,s)%Le noeud inhibiteur est présent
                                Pinh=Pinh*(1-inh(1,Gi(j,4)));%Probabilité que l'inhib j>i échoue
                            end
                        end
                    end
                    %Application de la proba (selon si apparition ou pas)
                    P=(Pini+(1-Pini)*(1-Pimp))*Pinh*Pact;
                else%Même procédure pour la survie
                    Pini=ini(2,I(i));%Présence à l'arrivée
                    Pimp=1;
                    Pinh=1;
                    Pact=act(2,A(i,t,s));
                    for j=1:size(Gi,1)
                        if Gi(j,3)==1
                            if D(Gi(j,1),t,s)
                                Pimp=Pimp*(1-imp(2,Gi(j,4)));%Probabilité que l'impulseur j>i échoue
                            end
                        end
                        if Gi(j,3)==2
                            if D(Gi(j,1),t,s)
                                Pinh=Pinh*(1-inh(2,Gi(j,4)));%Probabilité que l'inhib j>i échoue
                            end
                        end
                    end
                    P=(Pini+(1-Pini)*(1-Pimp))*Pinh*Pact;
                end
                D(i,t+1,s)=rand()<P;
            end
        end
    end