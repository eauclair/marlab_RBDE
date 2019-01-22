%fonction qui calcule la partie de la loglike du SBM qui depend de alpha
%% INPUT
%G : matrice d'adjacence du graphe, N*N*3
% TL : le vecteur colonne des niveaux trophiques des esp?ces
% alpha : param?tre du SBM qui modelise les interactions +
%% OUTPUT
% Y : les termes dans logP(LG | psi) qui dependent de alpha

%%
function [ Y ] = lvSBMalpha( Gimp, TL,  alpha)


n = size(Gimp,1); %nb de sommets
%TLp = TL';
A = repmat(TL,1,n);
B = A';
delta = A -B;

test = ones(n,n);
test(delta >-1) = 0;
testft = ones(n,n);
testft(find(TL==0),:) = 0;
testft(:,find(TL==0)) = 0;
test=test.*testft;

Y = alpha*test.*delta.*Gimp;

lo = exp(alpha*delta);
lo = log(1 + lo);

 Y = Y - lo.*test;
 


Y = sum(sum(Y));

end

