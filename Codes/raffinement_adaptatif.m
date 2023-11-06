function [Nbpt2,Nbtri2,Coorneu2,Refneu2,Numtri2,s,est] = raffinement_adaptatif(Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Na,Barycentres,Numaretes,Nbaretes_int,Tri_ar,Ar_tri,Numaretes_bords,homogene)

if strcmp(homogene,'oui')
    [U_h,Nbtri,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,KK] = test(Ar_tri,Nbaretes_int,Na,Numaretes,Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Numaretes_bords);
    G_U(:,1) = pi*cos(pi*Barycentres(:,1)).*sin(pi*Barycentres(:,2));
    G_U(:,2) = pi*sin(pi*Barycentres(:,1)).*cos(pi*Barycentres(:,2));
else
    [U_h,Nbtri,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,KK] = test_non_h(Barycentres,Tri_ar,Ar_tri,Nbaretes_int,Na,Numaretes,Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Numaretes_bords);
    G_U = zeros(Na,2);
    for i=1:Na
        G_U(i,:) = grad_polaire(Barycentres(i,1),Barycentres(i,2))';
    end
end

Aire_tri = Aire_triangles(Coorneu,Nbtri,Numtri);

estimateurs = estimation(U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,Barycentres,Aire_tri,homogene);
est = sqrt(sum(estimateurs));


s = 0;
erreurs = zeros(Nbtri,1);
for K=1:Nbtri
    a = 0;
    for F=Tri_ar(K,:)
        vec = zeros(2,1);
        vec(1) = (U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),2)-Coorneu(Numtri(K,1),2)) + U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),2)-Coorneu(Numtri(K,2),2)) + U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),2)-Coorneu(Numtri(K,3),2)))/Aire_tri(K);
        vec(2) = (-U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),1)-Coorneu(Numtri(K,1),1)) - U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),1)-Coorneu(Numtri(K,2),1)) - U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),1)-Coorneu(Numtri(K,3),1)))/Aire_tri(K);
        a = a + (G_U(F,:)-vec')*(G_U(F,:)'-vec);
    end
    erreurs(K) = a*Aire_tri(K)/3;
    s = s + a*Aire_tri(K)/3;
end


%{
figure;
for K=1:Nbtri
    S1 = Numtri(K,1);
    S2 = Numtri(K,2);
    S3 = Numtri(K,3);
    T = [1 2 3];
    X = [Coorneu(S1,1);Coorneu(S2,1);Coorneu(S3,1)];
    Y = [Coorneu(S1,2);Coorneu(S2,2);Coorneu(S3,2)];
    Z = [sqrt(estimateurs(K));sqrt(estimateurs(K));sqrt(estimateurs(K))];
    trisurf(T,X,Y,Z);
    hold on;
end
view(2);
colorbar;

figure;
for K=1:Nbtri
    S1 = Numtri(K,1);
    S2 = Numtri(K,2);
    S3 = Numtri(K,3);
    T = [1 2 3];
    X = [Coorneu(S1,1);Coorneu(S2,1);Coorneu(S3,1)];
    Y = [Coorneu(S1,2);Coorneu(S2,2);Coorneu(S3,2)];
    Z = [sqrt(erreurs(K));sqrt(erreurs(K));sqrt(erreurs(K))];
    trisurf(T,X,Y,Z);
    hold on;
end
view(2);
colorbar;
%}

estimateurs = estimateurs/sum(estimateurs);
classement = [estimateurs (1:Nbtri)'];
classement = sortrows(classement,1,'descend');
er = 0;
compteur = 1;
triangles_a_traiter = zeros(Nbtri,1);
while compteur<=Nbtri && er<0.5
    er = er + classement(compteur,1);
    triangles_a_traiter(compteur) = classement(compteur,2);
    compteur = compteur + 1;
end
indicatrice = zeros(Nbtri,1);
voisins = zeros(Nbtri,4);
faces_a_traiter = zeros(Na,1);
for K=triangles_a_traiter'
    if K~=0
        indicatrice(K) = 1;
    end
end
for K=1:Nbtri
    if indicatrice(K) == 0
        for F=Tri_ar(K,:)
            faces_a_traiter(F) = 1;
            if Ar_tri(F,1) ~= K
                K2 = Ar_tri(F,1);
            elseif Ar_tri(F,2) ~= 0
                K2 = Ar_tri(F,2);
            else
                K2 = 0;
            end
            if K2~=0 && indicatrice(K2)==1
                voisins(K,1) = voisins(K,1) + 1;
                voisins(K,voisins(K,1)+1) = K2;
            end
        end
    else
        voisins(K,1) = -1;
    end
end
l0 = find(voisins(:,1)==0);
l1 = find(voisins(:,1)==1);
l2 = find(voisins(:,1)==2);
l3 = find(voisins(:,1)==3);
l_traiter = find(indicatrice==1);
Nb_voisins0 = length(l0);
Nb_voisins1 = length(l1);
Nb_voisins2 = length(l2);
Nb_voisins3 = length(l3);
Nb_a_traiter = length(l_traiter);
Nbtri2 = 4*(Nb_a_traiter + Nb_voisins3) + 3*Nb_voisins2 + 2*Nb_voisins1 + Nb_voisins0;
Numtri2 = zeros(Nbtri2,3);
compteur = 1;
compteur_aretes = 1;
nv_aretes = zeros(Na,1);
for K=[l_traiter' l3']
    S1 = Numtri(K,1);
    S2 = Numtri(K,2);
    S3 = Numtri(K,3);
    Ar1 = Tri_ar(K,1);
    Ar2 = Tri_ar(K,2);
    Ar3 = Tri_ar(K,3);
    bar = [0 0 0];
    if nv_aretes(Ar1)==0
        bar(1) = Nbpt + compteur_aretes;
        nv_aretes(Ar1) = Nbpt + compteur_aretes;
        compteur_aretes = compteur_aretes + 1;
    else
        bar(1) = nv_aretes(Ar1);
    end
    if nv_aretes(Ar2)==0
        bar(2) = Nbpt + compteur_aretes;
        nv_aretes(Ar2) = Nbpt + compteur_aretes;
        compteur_aretes = compteur_aretes + 1;
    else
        bar(2) = nv_aretes(Ar2);
    end
    if nv_aretes(Ar3)==0
        bar(3) = Nbpt + compteur_aretes;
        nv_aretes(Ar3) = Nbpt + compteur_aretes;
        compteur_aretes = compteur_aretes + 1;
    else
        bar(3) = nv_aretes(Ar3);
    end
    Numtri2(4*compteur-3,:) = [S1 bar(1) bar(3)];
    Numtri2(4*compteur-2,:) = [bar(1) S2 bar(2)];
    Numtri2(4*compteur-1,:) = [bar(2) S3 bar(3)];
    Numtri2(4*compteur,:) = [bar(1) bar(2) bar(3)];
    compteur = compteur + 1;
end
compteur = 4*(Nb_a_traiter + Nb_voisins3) + 1;
for K=l2'
    Ar1 = Tri_ar(K,1);
    Ar2 = Tri_ar(K,2);
    Ar3 = Tri_ar(K,3);
    for F=Tri_ar(voisins(K,2),:)
        if F==Ar1 || F==Ar2 || F==Ar3
            arete1 = F;
        end
    end
    for F=Tri_ar(voisins(K,3),:)
        if F==Ar1 || F==Ar2 || F==Ar3
            arete2 = F;
        end
    end
    for i=1:3
        if Tri_ar(K,i) == arete1
            indice1 = i;
        end
        if Tri_ar(K,i) == arete2
            indice2 = i;
        end
    end
    if mod(indice1-indice2,3)==1
        temp = indice1;
        indice1 = indice2;
        indice2 = temp;
    end %indice1+1=indice2 modulo3
    Numtri2(compteur,:) = [nv_aretes(Tri_ar(K,indice2)) Numtri(K,mod(indice2,3)+1) Numtri(K,mod(indice1+2,3)+1)];
    compteur = compteur + 1;
    Numtri2(compteur,:) = [nv_aretes(Tri_ar(K,indice2)) Numtri(K,mod(indice2+1,3)+1) nv_aretes(Tri_ar(K,indice1))];
    compteur = compteur + 1;
    Numtri2(compteur,:) = [nv_aretes(Tri_ar(K,indice2)) nv_aretes(Tri_ar(K,indice1)) Numtri(K,indice2)];
    compteur = compteur + 1;
end
for K=l1'
    Ar1 = Tri_ar(K,1);
    Ar2 = Tri_ar(K,2);
    Ar3 = Tri_ar(K,3);
    for F=Tri_ar(voisins(K,2),:)
        if F==Ar1 || F==Ar2 || F==Ar3
            arete = F;
        end
    end
    for i=1:3
        if Numtri(K,i)~=Numaretes(arete,1) && Numtri(K,i)~=Numaretes(arete,2)
            indice_sommet = i;
        end
    end
    Numtri2(compteur,:) = [Numtri(K,indice_sommet) Numtri(K,mod(indice_sommet,3)+1) nv_aretes(arete)];
    compteur = compteur + 1;
    Numtri2(compteur,:) = [Numtri(K,indice_sommet) nv_aretes(arete) Numtri(K,mod(indice_sommet+1,3)+1)];
    compteur = compteur + 1;
end
for K=l0'
    Numtri2(compteur,:) = Numtri(K,:);
    compteur = compteur + 1;
end
Nbpt2 = Nbpt + compteur_aretes - 1;
Coorneu2 = zeros(Nbpt2,2);
Coorneu2(1:Nbpt,:) = Coorneu(1:Nbpt,:);
Refneu2 = zeros(Nbpt2,1);
Refneu2(1:Nbpt) = Refneu(1:Nbpt);
for F=1:Na
    if nv_aretes(F)~=0
        Coorneu2(nv_aretes(F),:) = (Coorneu(Numaretes(F,1),:) + Coorneu(Numaretes(F,2),:))/2;
        if Numaretes(F,3) == 2
            Refneu2(nv_aretes(F)) = 1;
        end
    end
end
%{
figure;
% maillage
trimesh(Numtri2(:,1:3),Coorneu2(:,1),Coorneu2(:,2),zeros(Nbpt2,1));
hold on;
plot(Coorneu2(:,1),Coorneu2(:,2),'.r', "markersize", 10);
view(2);
%}
