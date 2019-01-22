

%% INPUT
% TL : vecteur colonne des niveaux trophiques des especes
%G : matrice d'adjacence du graphe n*n*3
%
%% OUTPUT
% estimation du parametre beta2 du modele SBM
%
% NP 23/01/17 correction bug formule

function [ beta2 ] = estimbeta2( Gimp, Ginh, TL )

n = size(Gimp,1); %nb de sommets
A = repmat(TL,1,n);
B = A';
delta = A -B;

test = ones(n,n);
test(delta > 0) = 0;

testft = ones(n,n);
testft(find(TL==0),:) = 0;
testft(:,find(TL==0)) = 0;
test=test.*testft;

if sum(sum(test.*(1-Gimp)))
    beta2 = (sum(sum(test.*Ginh.*(1-Gimp)))) / sum(sum(test.*(1-Gimp)));
else
    beta2=0;
end
end
