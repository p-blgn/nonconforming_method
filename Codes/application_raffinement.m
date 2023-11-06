homogene = 'oui'
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes]=lecture_msh("geomCarre.msh");
nombre_de_raffinements = 12;
S = zeros(nombre_de_raffinements,1);
E = zeros(nombre_de_raffinements,1);
Dim = zeros(nombre_de_raffinements,1);
Nbaretes_ext = 0;%n'est pas utilisé
Numaretes_bords = 0;%idem
regularite = zeros(nombre_de_raffinements,1);
for i=1:nombre_de_raffinements
    i
    [Barycentres,Na,Numaretes,isavertex,Ar_tri,Nbaretes_int,Tri_ar] = Initialisation_matrices(Coorneu,Nbtri,Numtri,Nbaretes_ext,Numaretes_bords,Nbpt);
    reg_max = 0;
    l_ar = longueurs(Na,Coorneu,Numaretes);
    Aire_tri = Aire_triangles(Coorneu,Nbtri,Numtri);
    for K=1:Nbtri
        rho = Aire_tri(K)/(l_ar(Tri_ar(K,1)) + l_ar(Tri_ar(K,2)) + l_ar(Tri_ar(K,3)));
        h = max([l_ar(Tri_ar(K,1)) l_ar(Tri_ar(K,2)) l_ar(Tri_ar(K,3))]);
        reg = h/rho;
        if reg>reg_max
            reg_max = reg;
        end
    end
    regularite(i) = reg_max;
    Dim(i) = Na;
    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,s,est] = raffinement_adaptatif(Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Na,Barycentres,Numaretes,Nbaretes_int,Tri_ar,Ar_tri,Numaretes_bords,homogene);
    S(i) = sqrt(s);
    E(i) = est;
end
figure;
loglog(Dim,S,'r-*');
hold on;
loglog(Dim,E,'b-*');
legend("Erreur","Estimateur");
figure;
semilogx(Dim,E./S,'b-*');
legend("Indice d'efficacité");
figure;
semilogx(Dim,regularite,'r-*');
legend("Régularité");