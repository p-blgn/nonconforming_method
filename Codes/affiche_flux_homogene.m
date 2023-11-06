function affiche_flux_homogene(nom_maillage)
[U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,Barycentres,Aire_tri,Ar_tri] = TP1(nom_maillage);
[grad,flux_points,points,centres] = reconstruction_flux_homogene(Numtri,Nbtri,Na,Barycentres,Tri_ar,U_h,Coorneu,Aire_tri,Ar_tri);
figure;
Nbpt = length(Coorneu);
trimesh(Numtri(:,1:3),Coorneu(:,1),Coorneu(:,2),zeros(Nbpt,1));
hold on;
quiver(points(:,1),points(:,2),flux_points(:,1),flux_points(:,2));
hold on;
quiver(centres(:,1),centres(:,2),-grad(:,1),-grad(:,2),"r");
legend("Maillage","\sigma_h","-\nabla.u_h");
view(2);
end