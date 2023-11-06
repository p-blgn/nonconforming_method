function [U_h_plein,Nbtri,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,KK] = test(Ar_tri,Nbaretes_int,Na,Numaretes,Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Numaretes_bords)

% declarations
% ------------

Aire_tri = Aire_triangles(Coorneu,Nbtri,Numtri);
Longueurs_Ar = longueurs(Na,Coorneu,Numaretes);
KK = sparse(Nbaretes_int, Nbaretes_int);
B = zeros(Nbaretes_int,1);
Ns = Nbpt;
compteur1 = 0;
for i=1:Na
    if Numaretes(i,3) == 1
        compteur1 = compteur1 + 1;
        compteur2 = 0;
       for j=1:Na
           if Numaretes(j,3) == 1 
               compteur2 = compteur2 + 1;
               if compteur1==compteur2
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
compteur = 1;
for F=1:Na
    if Numaretes(F,3) == 1
        s = 0;
        for K=1:2
            if Ar_tri(F,K)~=0
                x = (Coorneu(Numaretes(F,1),1) + Coorneu(Numaretes(F,2),1))/2;
                y = (Coorneu(Numaretes(F,1),2) + Coorneu(Numaretes(F,2),2))/2;
                s = s + Aire_tri(Ar_tri(F,K))*f(x,y)/3;
            end
        end
        B(compteur) = s;
        compteur = compteur + 1;
    end
end
U_h = KK\B;
U_h_plein = zeros(Na,1);
compteur = 0;
for F=1:Na
    if Numaretes(F,3)==1
        compteur = compteur + 1;
        U_h_plein(F) = U_h(Numaretes(F,4));
    else
        U_h_plein(F) = 0;
    end
end
end