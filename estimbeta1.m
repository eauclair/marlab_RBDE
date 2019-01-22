%% INPUT
% TL : vecteur colonne des niveaux trophiques des especes
%G : matrice d'adjacence du graphe n*n*3
%
%% OUTPUT
% estimation du parametre beta1 du modele SBM
%

function [ beta1 ] = estimbeta1( Ginh, TL )

n = size(Ginh,1); %nb de sommets
A = repmat(TL,1,n);
B = A';
delta = A -B;

test = ones(n,n);
test(delta < 1) = 0;

testft = ones(n,n);
testft(find(TL==0),:) = 0;
testft(:,find(TL==0)) = 0;
test=test.*testft;

if sum(sum(test))
    beta1 = (sum(sum(test.*Ginh))) / sum(sum(test));
else
    beta1=0;
end
end

