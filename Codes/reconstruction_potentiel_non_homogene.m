function Val_sommets2 = reconstruction_potentiel_non_homogene(Val_sommets,Numtri,Ns,Nbtri,Refneu,Coorneu)
maxi = 0;
Nb_triangles_par_sommet = zeros(Ns,1);
Val_sommets2 = zeros(Ns,1);
for K=1:Nbtri
    for i=1:3
        S = Numtri(K,i);
        Nb_triangles_par_sommet(S) = Nb_triangles_par_sommet(S) + 1;
    end
end
maxi = max(Nb_triangles_par_sommet);
Som_tri = zeros(Ns, maxi + 1);
for K=1:Nbtri
    for i=1:3
        S = Numtri(K,i);
        Som_tri(S,1) = Som_tri(S,1) + 1; %Sert à trouver la colonne pour mettre le nouveau triangle
        Som_tri(S,Som_tri(S,1) + 1) = Val_sommets(K,i);
    end
end
for S=1:Ns
    if Refneu(S)==0
        Val_sommets2(S) = sum(Som_tri(S,2:Som_tri(S,1)+1))/Som_tri(S,1);
    else
        Val_sommets2(S) = u_bord(Coorneu(S,1),Coorneu(S,2));
    end
end