%% INPUT
% TL : vecteur colonne des niveaux trophiques des especes
%G : matrice d'adjacence du graphe n*n*3
%
%% OUTPUT
% estimation du parametre beta1 du modele SBM
%

function [ alphaft betaft] = estimft( Gimp, Ginh, TL )

n = size(Gimp,1); %nb de sommets

test = zeros(n,n);
test(find(TL==0),:) = 1;
test(:,find(TL==0)) = 1;

nFt=sum(sum(test));
if nFt
    alphaft = (sum(sum(test.*Gimp))) / nFt;
    betaft = (sum(sum(test.*Ginh))) / nFt;
else
    alphaft=0.01;
    betaft=0.01;
end
    
end

