%ltoA(L)
%à partir d'un tableau d'arcs A(espèce1, espece2, relation), renvoie une matrice d'adjacence
%rempli selon le type de relation
function[Y]=ltoA(L,n,nbE)
    Y=zeros(n,n,nbE);
    for i=1:size(L,1)
        Y(L(i,1),L(i,2),L(i,3))=1;
    end