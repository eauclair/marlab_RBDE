%Calcul de la log-vraisemblance d'une transition particulière
%d'un noeud i de t vers t+1 dans le jeu de données s
function[P]=lvTr(D, Gi, A, i, t, s, iniI,imp,inh,act)
%Absence au départ
if D(i,t,s)==0
    Pini=iniI(1);%Apparition spontanée (On ne donne en entrée que les paramètres d'initialisation relatifs au noeud i)
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
                %Finh=1;
            end
        end
    end
    %Application de la proba (selon si apparition ou pas)
    if D(i,t+1,s)
        P=(Pini+(1-Pini)*(1-Pimp))*Pinh*Pact;
    else
        P=1-(Pini+(1-Pini)*(1-Pimp))*Pinh*Pact;
    end
else%Même procédure pour la survie
    Pini=iniI(2);%Présence à l'arrivée
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
    if D(i,t+1,s)
        P=(Pini+(1-Pini)*(1-Pimp))*Pinh*Pact;
    else
        P=1-(Pini+(1-Pini)*(1-Pimp))*Pinh*Pact;
    end
end
P=log(P);