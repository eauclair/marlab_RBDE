function[G]=GSimulSBM(XSBM, TL, k)
n=size(TL, 2);

% la partie de la logvraisemblance qui vient du SBM
alpha = XSBM(1);
beta1 = XSBM(2);
beta2 = XSBM(3);
alphaft = XSBM(4);
betaft = XSBM(5);
G=zeros(0,4);
deg=zeros(n,1);%degré de chaque noeud
for i=1:n
    for j = 1 :n
        if deg(j)<k
            if j~=i
                deltaji = TL(j) - TL(i);
                isft=TL(j)*TL(i);
                isft=~isft;
                if rand()<(alpha*deltaji)*(deltaji >0) + (alphaft)*isft
                    G=[G;i j 1 1];
                    deg(j)=deg(j)+1;
                else
                    if rand()<(beta2)*(deltaji >=0) +  (beta1)*(deltaji<0) + (betaft)*isft
                        G=[G;i j 2 1];
                        deg(j)=deg(j)+1;
                    end
                end
                %C(Gi(j,3)) = (log(beta2) -log(1-beta2))*(deltaji <=0) +  (log(beta1) - log(1-beta1))*(deltaji>0) + (log(betaft) -log(1-betaft))*isft;
                %C(Gfn(j))  = -(log(beta2) -log(1-beta2))*(deltaji <0);
            end
        end
    end
end