function estimateurs = estimation(U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,Barycentres,Aire_tri,homogene)

if strcmp(homogene,'oui')
    Val_sommets = vers_sommets(U_h,Nbtri,Tri_ar,Na,Numaretes);
    Val_sommets2 = reconstruction_potentiel(Val_sommets,Numtri,Ns,Nbtri,Refneu);
else
    Val_sommets = vers_sommets_non_homogene(U_h,Nbtri,Tri_ar,Na,Numaretes);
    Val_sommets2 = reconstruction_potentiel_non_homogene(Val_sommets,Numtri,Ns,Nbtri,Refneu,Coorneu);
end

estimateurs = zeros(Nbtri,1);
Val_bar = zeros(Na,1); %valeur de s_h sur les barycentres, on peut alors appliquer la formule pour avoir le gradient
for F=1:Na
    Val_bar(F) = (Val_sommets2(Numaretes(F,1)) + Val_sommets2(Numaretes(F,2)))/2;
end

for K=1:Nbtri
    B = [Barycentres(Tri_ar(K,1),:); Barycentres(Tri_ar(K,2),:); Barycentres(Tri_ar(K,3),:)];
    centre = (B(1,:) + B(2,:) + B(3,:))/3;
    moyenne = (f(B(1,1),B(1,2)) + f(B(2,1),B(2,2)) + f(B(3,1),B(3,2)))/3;
    if strcmp(homogene,'oui')
        e1 = (moyenne^2)/12*Aire_tri(K)*(dot(B(1,:)-centre,B(1,:)-centre) + dot(B(2,:)-centre,B(2,:)-centre) + dot(B(3,:)-centre,B(3,:)-centre));%cas homogène
    else
        e1 = 0;
    end
    grad_u(1) = (U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),2)-Coorneu(Numtri(K,1),2)) + U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),2)-Coorneu(Numtri(K,2),2)) + U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),2)-Coorneu(Numtri(K,3),2)))/Aire_tri(K);
    grad_u(2) = (-U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),1)-Coorneu(Numtri(K,1),1)) - U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),1)-Coorneu(Numtri(K,2),1)) - U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),1)-Coorneu(Numtri(K,3),1)))/Aire_tri(K);
    grad_s(1) = (Val_bar(Tri_ar(K,1))*(Coorneu(Numtri(K,2),2)-Coorneu(Numtri(K,1),2)) + Val_bar(Tri_ar(K,2))*(Coorneu(Numtri(K,3),2)-Coorneu(Numtri(K,2),2)) + Val_bar(Tri_ar(K,3))*(Coorneu(Numtri(K,1),2)-Coorneu(Numtri(K,3),2)))/Aire_tri(K);
    grad_s(2) = (-Val_bar(Tri_ar(K,1))*(Coorneu(Numtri(K,2),1)-Coorneu(Numtri(K,1),1)) - Val_bar(Tri_ar(K,2))*(Coorneu(Numtri(K,3),1)-Coorneu(Numtri(K,2),1)) - Val_bar(Tri_ar(K,3))*(Coorneu(Numtri(K,1),1)-Coorneu(Numtri(K,3),1)))/Aire_tri(K);
    e2 = Aire_tri(K)*dot(grad_u-grad_s,grad_u-grad_s);
    estimateurs(K) = e1 + e2;
end
end