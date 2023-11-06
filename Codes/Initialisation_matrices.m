function [Barycentres,Na,Numaretes,isavertex,Ar_tri,Nbaretes_int,Tri_ar] = Initialisation_matrices(Coorneu,Nbtri,Numtri,Nbaretes_ext,Numaretes_bords,Nbpt)

isavertex = zeros(Nbpt, 1); % tableau booleen "is a vertex" reperant les noeuds qui sont des sommets (et pas des milieux)
M_ar = sparse(Nbpt,Nbpt);
Numaretes = [];
Barycentres = [];
Ar_tri = [];
Na = 0;
Nbaretes_int = 0;
Tri_ar = zeros(Nbtri, 3);
%On traite les arêtes de chaque triangle
for l=1:Nbtri
    w = 1;
    for i=1:3 %[1,2] [2,3] [3,1] pour garder le même sens pour les arêtes dans le triangle
        j = mod(i,3)+1;
            if M_ar(Numtri(l,i),Numtri(l,j))==0 %Si l'arête n'a pas encore été traitée
                Na = Na + 1;  %On augmente de nombre d'arêtes
                Ar_tri = [Ar_tri; 0 0];
                M_ar(Numtri(l,i),Numtri(l,j)) = Na; %On remplit la matrice de manière symétrique
                M_ar(Numtri(l,j),Numtri(l,i)) = Na;
                Numaretes = [Numaretes; Numtri(l,i) Numtri(l,j) 2 0]; %Numaretes contient les deux sommets extrémités, 2 si l'arête est sur le bord et 1 sinon et la dernière colonne contient sa numérotation dans le système int/ext
                Barycentres = [Barycentres;(Coorneu(Numaretes(Na,1),1) + Coorneu(Numaretes(Na,2),1))/2 (Coorneu(Numaretes(Na,1),2) + Coorneu(Numaretes(Na,2),2))/2];
            end
            if M_ar(Numtri(l,i),Numtri(l,j))~=0
                if Ar_tri(M_ar(Numtri(l,i),Numtri(l,j)),1) == 0
                    Ar_tri(M_ar(Numtri(l,i),Numtri(l,j)),1) = l;
                else
                    Ar_tri(M_ar(Numtri(l,i),Numtri(l,j)),2) = l;
                    Numaretes(M_ar(Numtri(l,i),Numtri(l,j)),3) = 1;
                    Nbaretes_int = Nbaretes_int + 1;
                end
                Tri_ar(l,w) = M_ar(Numtri(l,i),Numtri(l,j));
                w = w + 1;
            end

    end
    isavertex(Numtri(l,1)) = 1;
    isavertex(Numtri(l,2)) = 1;
    isavertex(Numtri(l,3)) = 1;
end % for l
compteur1 = 0;
compteur2 = 0;
for F=1:Na
    if Numaretes(F,3)==1
        compteur1 = compteur1 + 1;
        Numaretes(F,4) = compteur1;
    else
        compteur2 = compteur2 + 1;
        Numaretes(F,4) = compteur2;
    end
end