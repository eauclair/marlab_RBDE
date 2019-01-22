function[mConfL] = compGraph(G0, G1, n, nImp, nInh)
%Matrice de confusion
mConfL=zeros(1,5);
for imp1=1:nImp
    mConfL=[mConfL;1 imp1 0 0 0];
    mConfL=[mConfL;0 0 1 imp1 0];
    for imp2=1:nImp
        mConfL=[mConfL;1 imp1 1 imp2 0];
    end
    for inh2=1:nInh
        mConfL=[mConfL;1 imp1 2 inh2 0];
    end
end
for inh1=1:nInh
    mConfL=[mConfL;2 inh1 0 0 0];
    mConfL=[mConfL;0 0 2 inh1 0];
    for imp2=1:nImp
        mConfL=[mConfL;2 inh1 1 imp2 0];
    end
    for inh2=1:nInh
        mConfL=[mConfL;2 inh1 2 inh2 0];
    end
end
for i=1:1:n
    for j=1:n
        G0ij=G0(G0(:,1)==i,:);
        G0ij=G0ij(G0ij(:,2)==j,:);
        if size(G0ij,1)>1
            error(strcat('Multiple edges ',{' '}, num2str(i),'>',num2str(j),' on G0'))
        end
        if isempty(G0ij)
           e0=[0 0];
        else
           e0=[G0ij(3) G0ij(4)];
        end
        G1ij=G1(G1(:,1)==i,:);
        G1ij=G1ij(G1ij(:,2)==j,:);
        if size(G1ij,1)>1
            error(strcat('Multiple edges ',{' '}, num2str(i),'>',num2str(j),' on G1'))
        end
        if isempty(G1ij)
            e1=[0 0];
        else
            e1=[G1ij(3) G1ij(4)];
        end
        L=find(mConfL(:,1)==e0(1));
        L=intersect(L,find(mConfL(:,2)==e0(2)));
        L=intersect(L,find(mConfL(:,3)==e1(1)));
        L=intersect(L,find(mConfL(:,4)==e1(2)));
        mConfL(L,5)=mConfL(L,5)+1;
    end
end
