function Val_sommets2 = reconstruction_potentiel(Val_sommets,Numtri,Ns,Nbtri,Refneu)
maxi = 0;
Nb_triangles_par_sommet = zeros(Ns,1);
Val_sommets2 = zeros(Ns,1);
for K=1:Nbtri
    for i=1:3
        S = Numtri(K,i);
        Nb_triangles_par_sommet(S) = Nb_triangles_par_sommet(S) + 1;
        Val_sommets2(S) = Val_sommets2(S) + Val_sommets(K,i);
    end
end
Val_sommets2 = Val_sommets2./Nb_triangles_par_sommet;
for S=1:Ns
    if Refneu(S)==1
        Val_sommets2(S) = 0;
    end
end