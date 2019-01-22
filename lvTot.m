function[lv]=lvTot(D, G, A, I,ini,imp,inh,act)
n=size(D,1);
T=size(D,2);
Q=size(D,3);
lv=0;
for i=1:n
    Gi=G(find(G(:,2)==i),:);
    iniI=ini(:,I(i));
    for t=1:T-1
        for s=1:Q
            lv=lv+lvTr(D,Gi, A, i, t, s,iniI,imp,inh,act);
            %test(i,t,s)=lvTr(D,Gi, A, i, t, s,iniI,imp,inh,act);
        end
    end
end