%lvILP_SBM(D, A, i, block,ini, imp, inh, act, XSBM, nbV, indVar, Ri)
function[C]=lvILP_SBM1(D, i, C, block, XSBM, indVar)
    %Variables pour output
    Gimp=indVar{1,1};
    Ginh=indVar{1,2};
    I_ini=indVar{1,3};
  % la partie de la logvraisemblance qui vient du SBM
    n = size(D,1);
    alpha = XSBM(1);
    beta1 = XSBM(2);
    beta2 = XSBM(3);
    alphaft = XSBM(4);
    betaft = XSBM(5);
    
    for j = 1 :n
        deltaji = block(j) - block(i);
        isft=block(j)*block(i);
        isft=~isft;
        C(Gimp) = (alpha*deltaji - log(1 + exp(alpha*deltaji)))*(deltaji < 0) + (log(alphaft) -log(1-alphaft))*isft ;
        C(Ginh) = (log(beta2))*(deltaji <=0) +  (log(beta1))*(deltaji>0) + (log(betaft))*isft;
        %ancienne formule, erreur
        %C(Gi(j,2)) = (alpha*deltaji - log(1 + exp(alpha*deltaji)))*(deltaji <0);
        %C(Gi(j,3)) =  (log(beta1) - log(1-beta1)*(deltaji>0)  + (log(beta2) -log(1-beta2)))*(deltaji <=0)  ; 
    end
%   C=-CindVar;