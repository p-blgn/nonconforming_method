function affiche_exo1_homogene_potentiel_reconstruit(nom_maillage)
[U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,~,~,~] = TP1(nom_maillage);
Val_sommets = vers_sommets(U_h,Nbtri,Tri_ar,Na,Numaretes);
Val_sommets2 = reconstruction_potentiel(Val_sommets,Numtri,Ns,Nbtri,Refneu);
figure;
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),Val_sommets2);
view(3);
shading interp
% shading faceted
% shading flat
colorbar;