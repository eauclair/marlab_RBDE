%% Expérimentation sur modèle "complexe"
par=8;
n=30;
T=15;
Q=10;
k0=5;
nBlock=5;
%Nombre de paramètres du modèle
nIni=2;%Nombre d'initialisations différentes
nImp=2;%Nombre d'impulseurs différents
nInh=2;%Nombre d'inhibiteurs différents
nAct=2;%Nombre d'actions différentes
nEt=nImp+nInh;%Nombre total d'étiquettes
nSBM=5;
script_param_basique_3
script_simu%Simulation de graphe/données
%Gsel = [i j type(inh ou imp) force present/absent]
Gsel0=zeros(0,5);
Gsel20=Gr(rand(size(Gr,1),1)<0.2,:);
dataRBDE={{ini imp inh act} SBMblock Gr A I D Gsel20};
parsave(strcat('expe_RBDE_simu_data_',num2str(nExp),'.mat'),dataRBDE);
Gsel=Gsel0;
script_apprentissage%Apprentissage de structure
resAppt_noKn={graphFinTot paramFinTot timeERFinTot lvFinTot};
parsave(strcat('expe_RBDE_simu_noKn_',num2str(nExp),'.mat'),resAppt_noKn);
script_apprentissage_SBM%Apprentissage de structure avec prior SBM
resAppt_SBM={graphFinTot paramFinTot timeERFinTot lvFinTot paramSBMFinTot};
parsave(strcat('expe_RBDE_simu_SBM_',num2str(nExp),'.mat'),resAppt_SBM);
Gsel=Gsel20;
script_apprentissage%Apprentissage de structure avec arcs connus
resAppt_brutKn={graphFinTot paramFinTot timeERFinTot lvFinTot};
parsave(strcat('expe_RBDE_simu_brutKn_',num2str(nExp),'.mat'),resAppt_brutKn);