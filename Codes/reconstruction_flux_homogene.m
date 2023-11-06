function [grad,flux_points,points,centres] = reconstruction_flux_homogene(Numtri,Nbtri,Na,Barycentres,Tri_ar,U_h,Coorneu,Aire_tri,Ar_tri)
grad = zeros(Nbtri,2);
flux_points = zeros(Na,2);
points = zeros(3*Nbtri,2);
centres = zeros(Nbtri,2);
for K=1:Nbtri
    grad(K,1) = (U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),2)-Coorneu(Numtri(K,1),2)) + U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),2)-Coorneu(Numtri(K,2),2)) + U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),2)-Coorneu(Numtri(K,3),2)))/Aire_tri(K);
    grad(K,2) = (-U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),1)-Coorneu(Numtri(K,1),1)) - U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),1)-Coorneu(Numtri(K,2),1)) - U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),1)-Coorneu(Numtri(K,3),1)))/Aire_tri(K);
    B = [Barycentres(Tri_ar(K,1),:); Barycentres(Tri_ar(K,2),:); Barycentres(Tri_ar(K,3),:)];
    centre = (B(1,:) + B(2,:) + B(3,:))/3;
    centres(K,:) = centre;
    S1 = Numtri(K,1);
    S2 = Numtri(K,2);
    S3 = Numtri(K,3);
    points(3*K-2,:) = 0.6*Coorneu(S1,:) + 0.2*Coorneu(S2,:) + 0.2*Coorneu(S3,:);
    points(3*K-1,:) = 0.2*Coorneu(S1,:) + 0.6*Coorneu(S2,:) + 0.2*Coorneu(S3,:);
    points(3*K,:) = 0.2*Coorneu(S1,:) + 0.2*Coorneu(S2,:) + 0.6*Coorneu(S3,:);
end
for K=1:Nbtri
    B = [Barycentres(Tri_ar(K,1),:); Barycentres(Tri_ar(K,2),:); Barycentres(Tri_ar(K,3),:)];
    centre = centres(K);
    moyenne = (f(B(1,1),(B(1,2))) + f(B(2,1),(B(2,2))) + f(B(3,1),(B(3,2))))*Aire_tri(K)/3;
    flux_points(3*K-2,:) = -grad(K,:) + moyenne/2*(points(3*K-2,:)-centre);
    flux_points(3*K-1,:) = -grad(K,:) + moyenne/2*(points(3*K-1,:)-centre);
    flux_points(3*K,:) = -grad(K,:) + moyenne/2*(points(3*K,:)-centre);
end
end