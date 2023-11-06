function [U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,B,Barycentres,Aire_tri,Ar_tri] = principal_exo1_non_homogene(nom_maillage)
%nom_maillage = 'geom_non_homogene.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,~,Nbaretes_ext,Numaretes_bords] = lecture_msh(nom_maillage);


% declarations
% ------------

Aire_tri = Aire_triangles(Coorneu,Nbtri,Numtri);
[Barycentres,Na,Numaretes,isavertex,Ar_tri,Nbaretes_int,Tri_ar] = Initialisation_matrices(Coorneu,Nbtri,Numtri,Nbaretes_ext,Numaretes_bords,Nbpt);
Longueurs_Ar = longueurs(Na,Coorneu,Numaretes);
% find(isavertex~=0) donne la liste des numeros des sommets
numsommets = find(isavertex==1);
num_int = find(Numaretes(:,3)==1);
num_ext = find(Numaretes(:,3)==2);
Ns = length(numsommets);   % nombre de sommets
KK = sparse(Nbaretes_int, Nbaretes_int);
G = zeros(Nbaretes_ext,1);
compteur1 = 0;
for i=1:Na
    if Numaretes(i,3) == 1 
        compteur1 = compteur1 + 1;
        compteur2 = 0;
       for j=1:Na
           if Numaretes(j,3) == 1 
               compteur2 = compteur2 + 1;
               if compteur1==compteur2
                   sommet1 = numsommets(Numaretes(i,1));
                   sommet2 = numsommets(Numaretes(i,2));
                   triangle1 = Ar_tri(i,1);
                   triangle2 = Ar_tri(i,2);
                   if triangle2 == 0 %inutile puisqu'on regarde les arêtes intérieures
                       KK(compteur1,compteur2) = ((Longueurs_Ar(i))^2)*(1/Aire_tri(triangle1));
                   else
                       KK(compteur1,compteur2) = ((Longueurs_Ar(i))^2)*((1/Aire_tri(triangle1))+(1/Aire_tri(triangle2)));
                   end
               else
                   triangle11 = Ar_tri(i,1);
                   triangle12 = Ar_tri(i,2);
                   triangle21 = Ar_tri(j,1);
                   triangle22 = Ar_tri(j,2);
                   if triangle11==triangle21
                       triangle_commun = triangle11;
                   elseif triangle11==triangle22
                       triangle_commun = triangle11;
                   elseif triangle12==triangle21
                       triangle_commun = triangle12;
                   elseif (triangle12~=0 && triangle12==triangle22)
                       triangle_commun = triangle12;
                   else
                       triangle_commun = 0;
                   end
                   if triangle_commun ~= 0
                       if Numaretes(i,1) == Numaretes(j,1)
                           sommet_commun = Numaretes(i,1);
                       elseif Numaretes(i,1) == Numaretes(j,2)
                           sommet_commun = Numaretes(i,1);
                       else
                           sommet_commun = Numaretes(i,2);
                       end
                       sommets = [];
                       for k=1:3
                           if Numtri(triangle_commun,k)~=sommet_commun
                               sommets = [sommets Numtri(triangle_commun,k)];
                           end
                       end
                       KK(compteur1,compteur2) = -(Coorneu(sommet_commun,:)-Coorneu(sommets(1),:))*(Coorneu(sommet_commun,:)-Coorneu(sommets(2),:))'/Aire_tri(triangle_commun);
                   end
               end
           end
       end
    end
end
KK2 = sparse(Nbaretes_int,Nbaretes_ext);
for j=1:Nbaretes_ext
    F2 = num_ext(j);
    K = Ar_tri(F2,1);
    cotes = Tri_ar(K,:);
    for F=cotes
        if Numaretes(F,3)==1
            i = Numaretes(F,4);
            if Numaretes(F,1) == Numaretes(F2,1)
                sommet_commun = Numaretes(F,1);
            elseif Numaretes(F,1) == Numaretes(F2,2)
                sommet_commun = Numaretes(F,1);
            else
                sommet_commun = Numaretes(F,2);
            end
            sommets = [];
            for k=1:3
                if Numtri(K,k)~=sommet_commun
                    sommets = [sommets Numtri(K,k)];
                end
            end
            KK2(i,j) = -(Coorneu(sommet_commun,:)-Coorneu(sommets(1),:))*(Coorneu(sommet_commun,:)-Coorneu(sommets(2),:))'/Aire_tri(K);
        end
    end
end
for i=1:Nbaretes_ext
    F = num_ext(i);
    G(i) = u_bord(Barycentres(F,1),Barycentres(F,2));
end
B = -KK2*G;
U_h_int = KK\B;
U_h = zeros(Na,1);
for F=1:Na
    if Numaretes(F,3)==1
        U_h(F) = U_h_int(Numaretes(F,4));
    else
        U_h(F) = G(Numaretes(F,4));
    end
end