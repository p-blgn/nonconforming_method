function affiche_exo1_non_homogene_potentiel_reconstruit(nom_maillage)
[U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,~,~,~,~] = principal_exo1_non_homogene(nom_maillage);
Val_sommets = vers_sommets_non_homogene(U_h,Nbtri,Tri_ar,Na,Numaretes);
Val_sommets2 = reconstruction_potentiel_non_homogene(Val_sommets,Numtri,Ns,Nbtri,Refneu,Coorneu);
figure;
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),Val_sommets2);
view(3);
shading interp
% shading faceted
% shading flat
colorbar;