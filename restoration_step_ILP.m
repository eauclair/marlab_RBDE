%%%%%%%%%%%%%%%%%%%%%%%
%%%Learning Algorithm%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%Entries%%%
%D:Dataset (n*T*s)
%A:Covariates (n*T*s)
%k:maximum degree
%nodes:Liste de noeuds connus
%ini,imp,inh,act valeur des paramètres
function[Gi,Ini]=restoration_step_ILP(Ao, Bo, Co, indVar, i)
%Variables pour output
Gimp=indVar{1,1};
Ginh=indVar{1,2};
I_ini=indVar{1,3};
n=size(Gimp,2);
nImp=size(Gimp,1);
nInh=size(Ginh,1);
nIni=size(I_ini,2);
%ILP
X=cplexbilp(Co, Ao, Bo);
Gi=zeros(0,4);
indG=0;
for j=1:n
    for indImp=1:nImp
        if X(Gimp(indImp,j))
            indG=indG+1;
            Gi(indG,1)=j;
            Gi(indG,2)=i;
            Gi(indG,3)=1;
            Gi(indG,4)=indImp;
        end
    end
    for indInh=1:nInh
        if X(Ginh(indInh,j))
            indG=indG+1;
            Gi(indG,1)=j;
            Gi(indG,2)=i;
            Gi(indG,3)=2;
            Gi(indG,4)=indImp;
        end
        
    end
end
for indI=1:nIni
    if X(I_ini(indI))
        Ini=indI;
    end
end