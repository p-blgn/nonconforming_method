function [grad,flux_barycentres,centres] = reconstruction_flux_non_homogene(Numtri,Nbtri,Na,Barycentres,Tri_ar,U_h,Coorneu,Aire_tri,Ar_tri,B_vec,Numaretes)
grad = zeros(Nbtri,2);
flux_points = zeros(Na,2);
points = zeros(3*Nbtri,2);
centres = zeros(Nbtri,2);
for K=1:Nbtri
    grad(K,1) = (U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),2)-Coorneu(Numtri(K,1),2)) + U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),2)-Coorneu(Numtri(K,2),2)) + U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),2)-Coorneu(Numtri(K,3),2)))/Aire_tri(K);
    grad(K,2) = (-U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),1)-Coorneu(Numtri(K,1),1)) - U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),1)-Coorneu(Numtri(K,2),1)) - U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),1)-Coorneu(Numtri(K,3),1)))/Aire_tri(K);
    B = [Barycentres(Tri_ar(K,1),:) ;Barycentres(Tri_ar(K,2),:); Barycentres(Tri_ar(K,3),:)];
    centre = (B(1,:) + B(2,:) + B(3,:))/3;
    centres(K,:) = centre;
end
flux_barycentres = -grad;
end