%% Calcul des précision et rappel
%Demande : Gapp : graphe appris, Gr : graphe "réel", n : nombre de noeuds,  
%nImp, nInh : Nombre d'impulseurs et d'inhibiteurs
mConf=compGraph(Gapp,Gr, n, nImp, nInh);
sameType=mConf(mConf(:,1)==mConf(:,3),:);%Appris du même type (peu importe la force de l'étiquette)
sameTypeEdge=sameType(sameType(:,1)~=0,:);%Appris du même type (sauf pas d'arc)
sameTypeImp=sameType(sameType(:,1)==1,:);%Appris du même type (imp)
sameTypeInh=sameType(sameType(:,1)==2,:);%Appris du même type (inh)
equal=sameType(sameType(:,2)==sameType(:,4),:);%Correctement appris : total
equalNone=equal(equal(:,1)==0,:);%Correctement appris : pas d'arc
equalEdge=equal(equal(:,1)~=0,:);%Correctement appris : arc présent
equalImp=equal(equal(:,1)==1,:);%Correctement appris : immulseurs
equalInh=equal(equal(:,1)==2,:);%Correctement appris : inhibiteurs
G0None=mConf(mConf(:,1)==0,:);%Appris comme pas d'arête en G0
G0Edge=mConf(mConf(:,1)~=0,:);%Appris comme avec arête en G0
G0Imp=mConf(mConf(:,1)==1,:);%Appris comme imp en G0
G0Inh=mConf(mConf(:,1)==2,:);%Appris comme inh en G0
G1None=mConf(mConf(:,3)==0,:);%Appris comme pas d'arête en G1
G1Edge=mConf(mConf(:,3)~=0,:);%Appris comme avec arête en G1
G1Imp=mConf(mConf(:,3)==1,:);%Appris comme imp en G1
G1Inh=mConf(mConf(:,3)==2,:);%Appris comme inh en G1

precision=sum(equalEdge,1)/sum(G0Edge,1);%Précision (Tout arc)
recall=sum(equalEdge,1)/sum(G1Edge,1);%Rappel (Tout arc)
precisionNone=sum(equalNone,1)/sum(G0None,1);%Précision (Pas d'arc)
recallNone=sum(equalNone,1)/sum(G1None,1);%Rappel (Pas d'arc)
precisionImp=sum(equalImp,1)/sum(G0Imp,1);%Précision (Impulseurs)
recallImp=sum(equalImp,1)/sum(G1Imp,1);%Rappel (Impulseurs)
precisionInh=sum(equalInh,1)/sum(G0Inh,1);%Précision (Inhibiteurs)
recallInh=sum(equalInh,1)/sum(G1Inh,1);%Rappel (Inhibiteurs)