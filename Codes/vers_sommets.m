function Val_sommets = vers_sommets(U_h_complet,Nbtri,Tri_ar,Na,Numaretes)
%Les fonctions sont de degré 1 donc on peut calculer la valeur en chaque
%sommet pour chaque triangle
Val_sommets = zeros(Nbtri,3);
compteur = 0;
for K=1:Nbtri
    Val_sommets(K,1) = U_h_complet(Tri_ar(K,1)) - U_h_complet(Tri_ar(K,2)) + U_h_complet(Tri_ar(K,3));
    Val_sommets(K,2) = -U_h_complet(Tri_ar(K,3)) + U_h_complet(Tri_ar(K,1)) + U_h_complet(Tri_ar(K,2));
    Val_sommets(K,3) = U_h_complet(Tri_ar(K,2)) + U_h_complet(Tri_ar(K,3)) - U_h_complet(Tri_ar(K,1));
end
end