%% Script : résultats des expérmimentations
%Données en entrée
nImp=1;
nInh=1;
n=20;
%Initialisation des tables
tabResNoKnOV = cell2table(cell(15,5), 'VariableNames', {'Exp', 'app', 'sur', 'lv', 'time'});
tabResBrutKnOV = cell2table(cell(15,5), 'VariableNames', {'Exp', 'app', 'sur', 'lv', 'time'});
tabResSBMOV = cell2table(cell(15,5), 'VariableNames', {'Exp', 'app', 'sur', 'lv', 'time'});
tabResNoKnPar = cell2table(cell(15,13), 'VariableNames', {'Exp', 'RepsApp', 'RrhoApp', 'RtauApp', 'EepsApp', 'ErhoApp', 'EtauApp', 'RepsSur', 'RrhoSur', 'RtauSur', 'EepsSur', 'ErhoSur', 'EtauSur'});
tabResBrutKnPar = cell2table(cell(15,13), 'VariableNames', {'Exp', 'RepsApp', 'RrhoApp', 'RtauApp', 'EepsApp', 'ErhoApp', 'EtauApp', 'RepsSur', 'RrhoSur', 'RtauSur', 'EepsSur', 'ErhoSur', 'EtauSur'});
tabResSBMPar = cell2table(cell(15,19), 'VariableNames', {'Exp', 'RepsApp', 'RrhoApp', 'RtauApp', 'EepsApp', 'ErhoApp', 'EtauApp', 'RepsSur', 'RrhoSur', 'RtauSur', 'EepsSur', 'ErhoSur', 'EtauSur', 'RA1', 'EA1','RB1', 'EB1', 'RB2', 'EB2'});
tabResNoKnApp = cell2table(cell(15,5), 'VariableNames', {'Exp', 'preImp', 'recImp', 'preInh', 'recInh'});
tabResBrutKnApp = cell2table(cell(15,5), 'VariableNames', {'Exp', 'preImp', 'recImp', 'preInh', 'recInh'});
tabResSBMApp = cell2table(cell(15,5), 'VariableNames', {'Exp', 'preImp', 'recImp', 'preInh', 'recInh'});
%parcours des expés
for i=1:15
    switch mod(i,3)
    case 1
        isApp=0;
        isSur=1;
    case 2
        isApp=1;
        isSur=0;
    case 0
        isApp=1;
        isSur=1;
    end
    %Chargement des données de base
    try 
        load(strcat('expe_RBDE_simu_data_',num2str(i),'.mat'));
        Gr0=res{3};%Graphe réel
        X=res{1};%Paramètres réels
        XSBM=[0.8 0.3 0.2];%Paramètres SBM
        Gsel=res{7};%Arêtes connues
        GrSel=setdiff(Gr0,Gsel(:,1:4),'rows');%Graphe réel sans les arês connues
    catch
        tabResNoKnOV{i,:}={i, isApp, isSur, nan, nan};
        tabResBrutKnOV{i,:}={i, isApp, isSur, nan, nan};
        tabResSBMOV{i,:}={i, isApp, isSur, nan, nan};
        tabResNoKnPar{i,:}={i, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
        tabResBrutKnPar{i,:}={i, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
        tabResSBMPar{i,:}={i, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
        tabResNoKnApp{i,:}={i, nan, nan, nan, nan};
        tabResBrutKnApp{i,:}={i, nan, nan, nan, nan};
        tabResSBMApp{i,:}={i, nan, nan, nan, nan};
    end
    %%Expe 'Sans connaissances'
    try
        load(strcat('expe_RBDE_simu_noKn_',num2str(i),'.mat'));
        %{graphFinTot paramFinTot timeERFinTot lvFinTot}
        Gapp0=res{1};
        Xapp=res{2};
        time=res{3};
        lv=res{4};
        tabResNoKnOV{i,:}={i, isApp, isSur, lv, time};
        tabResNoKnPar{i,:}={i, X{1}(1), X{2}(1), X{3}(1), Xapp(1,1), Xapp(1,2), Xapp(1,3), X{1}(2), X{2}(2), X{3}(2), Xapp(2,1), Xapp(2,2), Xapp(2,3)};
        Gr=Gr0;
        script_preRec
        tabResNoKnApp{i,:}={i, precisionImp, recallImp, precisionInh, recallInh};
    catch
        tabResNoKnOV{i,:}={i, isApp, isSur, nan, nan};
        tabResNoKnPar{i,:}={i, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
        tabResNoKnApp{i,:}={i, nan, nan, nan, nan};
    end
    %%Expe 'Connaissances brutes'
    try
        load(strcat('expe_RBDE_simu_brutKn_',num2str(i),'.mat'));
        Gapp=res{1};
        Xapp=res{2};
        time=res{3};
        lv=res{4};
        tabResBrutKnOV{i,:}={i, isApp, isSur, lv, time};
        tabResBrutKnPar{i,:}={i, X{1}(1), X{2}(1), X{3}(1), Xapp(1,1), Xapp(1,2), Xapp(1,3), X{1}(2), X{2}(2), X{3}(2), Xapp(2,1), Xapp(2,2), Xapp(2,3)};
        Gr=GrSel;
        Gapp=setdiff(Gapp,Gsel(:,1:4),'rows');%Graphe appris sans les arês connues
        script_preRec
        tabResBrutKnApp{i,:}={i, precisionImp, recallImp, precisionInh, recallInh};
    catch
        tabResBrutKnOV{i,:}={i, isApp, isSur, nan, nan};
        tabResBrutKnPar{i,:}={i, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
        tabResBrutKnApp{i,:}={i, nan, nan, nan, nan};
    end
     %%Expe 'SBM'
    try
        load(strcat('expe_RBDE_simu_SBM_',num2str(i),'.mat'));
        Gapp=res{1};
        Xapp=res{2};
        time=res{3};
        lv=res{4};
        XSBMApp=res{5}(1:3);
        tabResSBMOV{i,:}={i, isApp, isSur, lv, time};
        tabResSBMPar{i,:}={i, X{1}(1), X{2}(1), X{3}(1), Xapp(1,1), Xapp(1,2), Xapp(1,3), X{1}(2), X{2}(2), X{3}(2), Xapp(2,1), Xapp(2,2), Xapp(2,3), XSBM(1), XSBMApp(1), XSBM(2), XSBMApp(2), XSBM(3), XSBMApp(3)};
        Gr=GrSel;
        Gapp=setdiff(Gapp,Gsel(:,1:4),'rows');%Graphe appris sans les arês connues
        script_preRec
        tabResSBMApp{i,:}={i, precisionImp, recallImp, precisionInh, recallInh};
    catch
        tabResSBMOV{i,:}={i, isApp, isSur, nan, nan};
        tabResSBMPar{i,:}={i, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
        tabResSBMApp{i,:}={i, nan, nan, nan, nan};
    end
end
writetable(tabResNoKnOV)
writetable(tabResNoKnOV)
writetable(tabResBrutKnOV)
writetable(tabResSBMOV)
writetable(tabResNoKnPar)
writetable(tabResBrutKnPar)
writetable(tabResSBMPar)
writetable(tabResNoKnApp)
writetable(tabResBrutKnApp)
writetable(tabResSBMApp)
