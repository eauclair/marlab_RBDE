%Estimation des paramètres : met les paramètres en un vecteur, les estime.
%Certains paramètres sont connus, d'autres non 
function[ini,imp,inh,act]=parEstKn(D, G, A, I,par0,parKn,parEq)
%Vecteur de tous les paramètres
Xapp=[];
Xsur=[];
%Vecteur des paramètres égaux
XappEq=[];
XsurEq=[];
%Vecteur des paramètres connus
XappKn=[];
XsurKn=[];
%Paramètres initiaux
ini0=par0{1};
imp0=par0{2};
inh0=par0{3};
act0=par0{4};
%Paramètres connus
iniKn=parKn{1};
impKn=parKn{2};
inhKn=parKn{3};
actKn=parKn{4};
%Paramètres égaux
iniEq=parEq{1};
impEq=parEq{2};
inhEq=parEq{3};
actEq=parEq{4};
%Nombre de paramètres de chaque type
nIni=size(ini0,2);
nImp=size(imp0,2);
nInh=size(inh0,2);
nAct=size(act0,2);
%Indices de chaque paramètre
%X=[iniapp1 ... iniappmax impapp1 .... impappmax inhapp actapp inisur impsur...]
lastApp=0;
lastSur=nIni+nImp+nInh+nAct;
last=0;
iniInd=zeros(2,nIni);
for i=1:nIni
    last=last+1;
    %Paramètres initiaux
    Xapp=[Xapp ini0(1,i)];
    Xsur=[Xsur ini0(2,i)];
    %Paramètres égaux
    XappEq=[XappEq iniEq(1,i)];
    XsurEq=[XsurEq iniEq(2,i)];
    %Paramètres connus
    XappKn=[XappKn iniKn(1,i)];
    XsurKn=[XsurKn iniKn(2,i)];
    %Indices de la grande matrice X
    iniInd(1,i)=lastApp+last;
    iniInd(2,i)=lastSur+last;
end
impInd=zeros(2,nImp);
for i=1:nImp
    last=last+1;
    %Paramètres initiaux
    Xapp=[Xapp imp0(1,i)];
    Xsur=[Xsur imp0(2,i)];
    %Paramètres égaux
    XappEq=[XappEq impEq(1,i)];
    XsurEq=[XsurEq impEq(2,i)];
    %Paramètres connus
    XappKn=[XappKn impKn(1,i)];
    XsurKn=[XsurKn impKn(2,i)];
    %Indices de la grande matrice X
    impInd(1,i)=lastApp+last;
    impInd(2,i)=lastSur+last; 
end
inhInd=zeros(2,nInh);
for i=1:nInh
    last=last+1;
    %Paramètres initiaux
    Xapp=[Xapp inh0(1,i)];
    Xsur=[Xsur inh0(2,i)];
    %Paramètres égaux
    XappEq=[XappEq inhEq(1,i)];
    XsurEq=[XsurEq inhEq(2,i)];
    %Paramètres connus
    XappKn=[XappKn inhKn(1,i)];
    XsurKn=[XsurKn inhKn(2,i)];
    %Indices de la grande matrice X
    inhInd(1,i)=lastApp+last;
    inhInd(2,i)=lastSur+last; 
end
actInd=zeros(2,nAct);
for i=1:nAct
    last=last+1;
    %Paramètres initiaux
    Xapp=[Xapp act0(1,i)];
    Xsur=[Xsur act0(2,i)];
    %Paramètres égaux
    XappEq=[XappEq actEq(1,i)];
    XsurEq=[XsurEq actEq(2,i)];
    %Paramètres connus
    XappKn=[XappKn actKn(1,i)];
    XsurKn=[XsurKn actKn(2,i)];
    %Indices de la grande matrice X
    actInd(1,i)=lastApp+last;
    actInd(2,i)=lastSur+last; 
end
%Vecteur des valeurs intiales
X0=[Xapp Xsur];
nPar=size(X0,2);
%Vecteur des paramètres égaux et connus
Xeq=[XappEq XsurEq];
Xkn=[XappKn XsurKn];
Aeq=[];
Beq=[];
for i=find(Xkn)
    Ai=zeros(1,nPar);
    Ai(i)=1;
    Aeq=[Aeq;Ai];
    Beq=[Beq X0(i)];
end
for i=unique(Xeq)
    grp=find(Xeq==i);
    ref=grp(1);
    grp=grp(2:size(grp,2));
    if size(grp>0)
        for j=grp
            Aj=zeros(1,nPar);
            Aj(ref)=1;
            Aj(j)=-1;
            Aeq=[Aeq;Aj];
            Beq=[Beq 0];
        end
    end
end
Beq=Beq';
%lvTot(D,G,A,I,ini0,imp0,inh0,act0)
f = @(X)-lvTot(D,G,A,I,reshape(X(iniInd),size(iniInd,1),size(iniInd,2)),reshape(X(impInd),size(impInd,1),size(impInd,2)),reshape(X(inhInd),size(inhInd,1),size(inhInd,2)),reshape(X(actInd),size(actInd,1),size(actInd,2)));
X=fmincon(f,X0,[],[],Aeq,Beq,zeros(size(X0)),ones(size(X0)),[],optimset('Display','off'));
ini=reshape(X(iniInd),size(iniInd,1),size(iniInd,2));
imp=reshape(X(impInd),size(impInd,1),size(impInd,2));
inh=reshape(X(inhInd),size(inhInd,1),size(inhInd,2));
act=reshape(X(actInd),size(actInd,1),size(actInd,2));