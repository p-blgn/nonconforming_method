function affiche_exo1_homogene(nom_maillage)
[U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,~,~,~] = TP1(nom_maillage);
Val_sommets = vers_sommets(U_h,Nbtri,Tri_ar,Na,Numaretes);
figure;
for K=1:Nbtri
    S1 = Numtri(K,1);
    S2 = Numtri(K,2);
    S3 = Numtri(K,3);
    T = [1 2 3];
    X = [Coorneu(S1,1);Coorneu(S2,1);Coorneu(S3,1)];
    Y = [Coorneu(S1,2);Coorneu(S2,2);Coorneu(S3,2)];
    Z = [Val_sommets(K,1);Val_sommets(K,2);Val_sommets(K,3)];
    trisurf(T,X,Y,Z);
    hold on;
end
view(3);
shading interp
% shading faceted
% shading flat
colorbar;