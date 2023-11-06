function affiche_flux_non_homogene(nom_maillage)
[U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,B_vec,Barycentres,Aire_tri,Ar_tri] = principal_exo1_non_homogene(nom_maillage);
[grad,flux_barycentres,centres] = reconstruction_flux_non_homogene(Numtri,Nbtri,Na,Barycentres,Tri_ar,U_h,Coorneu,Aire_tri,Ar_tri,B_vec,Numaretes);
figure;
Nbpt = length(Coorneu);
trimesh(Numtri(:,1:3),Coorneu(:,1),Coorneu(:,2),zeros(Nbpt,1));
hold on;
quiver(centres(:,1),centres(:,2),-grad(:,1),-grad(:,2),"r");
legend("Maillage","-\nabla.u_h = \sigma_h");
view(2);
end