function Val_sommets = vers_sommets_non_homogene(U_h,Nbtri,Tri_ar,Na,Numaretes);
%Les fonctions sont linéaires donc on peut calculer la valeur en chaque
%sommet pout chaque triangle
Val_sommets = zeros(Nbtri,3);
for K=1:Nbtri
    Val_sommets(K,1) = U_h(Tri_ar(K,1)) - U_h(Tri_ar(K,2)) + U_h(Tri_ar(K,3));
    Val_sommets(K,2) = -U_h(Tri_ar(K,3)) + U_h(Tri_ar(K,1)) + U_h(Tri_ar(K,2));
    Val_sommets(K,3) = U_h(Tri_ar(K,2)) + U_h(Tri_ar(K,3)) - U_h(Tri_ar(K,1));
end