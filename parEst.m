%Estimation des paramètres : met les paramètres en un vecteur, les estime, 
function[ini,imp,inh,act]=parEst(D, G, A, I,ini0,imp0,inh0,act0)
%Vecteur de tous les paramètres
Xapp=[];
Xsur=[];
%Indices de chaque paramètre
%X=[iniapp1 ... iniappmax impapp1 .... impappmax inhapp actapp inisur impsur...]
lastApp=0;
lastSur=size(ini0,2)+size(imp0,2)+size(inh0,2)+size(act0,2);
last=0;
for i=1:size(ini0,2)
    last=last+1;
    Xapp=[Xapp ini0(1,i)];
    Xsur=[Xsur ini0(2,i)];
    iniInd(1,i)=lastApp+last;
    iniInd(2,i)=lastSur+last; 
end
for i=1:size(imp0,2)
    last=last+1;
    Xapp=[Xapp imp0(1,i)];
    Xsur=[Xsur imp0(2,i)];
    impInd(1,i)=lastApp+last;
    impInd(2,i)=lastSur+last; 
end
for i=1:size(inh0,2)
    last=last+1;
    Xapp=[Xapp inh0(1,i)];
    Xsur=[Xsur inh0(2,i)];
    inhInd(1,i)=lastApp+last;
    inhInd(2,i)=lastSur+last; 
end
for i=1:size(act0,2)
    last=last+1;
    Xapp=[Xapp act0(1,i)];
    Xsur=[Xsur act0(2,i)];
    actInd(1,i)=lastApp+last;
    actInd(2,i)=lastSur+last; 
end
X0=[Xapp Xsur];
X0(:)=0.01;
%lvTot(D,G,A,I,ini0,imp0,inh0,act0)
f = @(X)-lvTot(D,G,A,I,reshape(X(iniInd),size(iniInd,1),size(iniInd,2)),reshape(X(impInd),size(impInd,1),size(impInd,2)),reshape(X(inhInd),size(inhInd,1),size(inhInd,2)),reshape(X(actInd),size(actInd,1),size(actInd,2)));
X=fmincon(f,X0,[],[],[],[],zeros(size(X0)),ones(size(X0)),[],optimset('Display','off'));
ini=reshape(X(iniInd),size(iniInd,1),size(iniInd,2));
imp=reshape(X(impInd),size(impInd,1),size(impInd,2));
inh=reshape(X(inhInd),size(inhInd,1),size(inhInd,2));
act=reshape(X(actInd),size(actInd,1),size(actInd,2));