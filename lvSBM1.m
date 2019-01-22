%fonction qui calcule la  loglike du SBM 
%% INPUT
% G : matrice d'adjacence du graphe N*N*3
% TL : le vecteur colonne des niveaux trophiques des especes
% XSBM  : [alpha, beta1, beta2, betaft, alphaft], les 5 params du SBM
%% OUTPUT
% Y : logP(LG | psi) 
function [ Y ] = lvSBM1(Gimp, Ginh, TL,  XSBM)


n = size(Gimp,1); %nb de sommets
A = repmat(TL,1,n);
B = A';
delta = A -B;
ft=find(TL==0);
%La valeur de TL est de 0.5 lorsqu'une des espèces est dans le bloc ft
delta(ft,:)=0.5;
delta(:,ft)=0.5;
 
%partie qui depend de alpha
alpha = XSBM(1);
test = ones(n,n);
test(delta >-1) = 0;

Y = alpha*test.*delta.*Gimp;

lo = exp(alpha*delta);
lo = log(1 + lo);

 Y = Y - lo.*test;
 
 %partie qui depend de beta1
 beta1 = XSBM(2);
 test = ones(n,n);
 test(delta < 1) = 0;
 
 toto  = log(beta1)*Ginh + log(1 - beta1)*(1- Ginh);
 toto = test.*toto;
 
 Y = Y + toto;
 
 %partie qui depend de beta2
 beta2 = XSBM(3);
 test = ones(n,n);
 test(delta > 0 ) = 0;
 
 toto  = log(beta2)*Ginh + log(1 - beta2)*(1- Ginh);
 toto = test.*toto;
 
 Y = Y + toto;
 
 if ~prod(TL)%On regarde s'il y a un bft
     %partie qui depend de alphaft et de betaft
     alphaft = XSBM(4);
     betaft = XSBM(5);
     test = ones(n,n);
     test(delta ~= 0.5) = 0;

     toto  = log(alphaft)*Gimp + log(1 - alphaft)*(1- Gimp);
     toto = test.*toto;
     Y = Y + toto;

     toto  = log(betaft)*Ginh + log(1 - betaft)*(1- Ginh);
     toto = test.*toto;
     Y = Y + toto;
 end
 
Y = sum(sum(Y));

end