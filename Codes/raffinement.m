function [Nbpt2,Nbtri2,Coorneu2,Refneu2,Numtri2,Nbaretes_ext2,Numaretes_bords2,Na2,Barycentres2,Numaretes2,Ar_tri2,Nbaretes_int2,Tri_ar2] = raffinement(Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Na,Barycentres,Numaretes,Nbaretes_int,Tri_ar)

Nbtri2 = 4*Nbtri;
Numtri2 = zeros(Nbtri2,3);
Nbpt2 = Na+Nbpt;
Coorneu2 = zeros(Nbpt2,2);
Coorneu2(1:Nbpt,:) = Coorneu;
Na2 = 2*Na+3*Nbtri;
Numaretes2 = zeros(Na2,4);
Refneu2 = zeros(Nbpt2,1);
Refneu2(1:Nbpt,:) = Refneu;
Nbaretes_ext2 = 2*Nbaretes_ext;
Nbaretes_int2 = Na2 - Nbaretes_ext2;
Ar_tri2 = zeros(Na2,2);
Tri_ar2 = zeros(Nbtri2,3);
Numaretes_bords2 = zeros(Nbaretes_ext2,2);
compteur = 0;
for F=1:Na
    Coorneu2(Nbpt+F,:) = Barycentres(F,:);
    Numaretes2(F,1) = Numaretes(F,1);
    Numaretes2(F,2) = Nbpt+F;
    Numaretes2(Na+F,1) = Numaretes(F,2);
    Numaretes2(Na+F,2) = Nbpt+F;
    Refneu2(Nbpt+F) = Numaretes(F,3)-1; %Numaretes(F,3) vaut 2 si arête du bord et 1 sinon
    Numaretes2(F,3) = Numaretes(F,3);
    Numaretes2(Na+F,3) = Numaretes(F,3);
    Numaretes2(F,4) = Numaretes(F,4);
    if Numaretes(F,3)==2
        compteur = compteur+2;
        Numaretes2(Na+F,4) = Numaretes(F,4) + Nbaretes_ext;
        Numaretes_bords2(compteur-1,:) = [Numaretes(F,1) Na+F];
        Numaretes_bords2(compteur,:) = [Na+F Numaretes(F,2)];
    else
        Numaretes2(Na+F,4) = Numaretes(F,4) + Nbaretes_int;
    end
end
for K=1:Nbtri
    S1 = Numtri(K,1);
    S2 = Numtri(K,2);
    S3 = Numtri(K,3);
    Ar1 = Tri_ar(K,1);
    Ar2 = Tri_ar(K,2);
    Ar3 = Tri_ar(K,3);
    bar = [0 0 0];
    bar(1) = Nbpt + Ar1;
    bar(2) = Nbpt + Ar2;
    bar(3) = Nbpt + Ar3;
    Numaretes2(2*Na+K*3-2,:) = [bar(1) bar(2) 1 2*Nbaretes_int+3*K-2];
    Numaretes2(2*Na+K*3-1,:) = [bar(2) bar(3) 1 2*Nbaretes_int+3*K-1];
    Numaretes2(2*Na+K*3,:) = [bar(1) bar(3) 1 2*Nbaretes_int+3*K];
    Numtri2(4*K-3,:) = [S1 bar(1) bar(3)];
    Numtri2(4*K-2,:) = [bar(1) S2 bar(2)];
    Numtri2(4*K-1,:) = [bar(2) S3 bar(3)];
    Numtri2(4*K,:) = [bar(1) bar(2) bar(3)];
    if Numaretes(Ar1,1)==S1
        Ar11 = Ar1;
        Ar12 = Ar1 + Na;
    else
        Ar12 = Ar1;
        Ar11 = Ar1 + Na;
    end
    if Numaretes(Ar2,1)==S2
        Ar21 = Ar2;
        Ar22 = Ar2 + Na;
    else
        Ar22 = Ar2;
        Ar21 = Ar2 + Na;
    end
    if Numaretes(Ar3,1)==S3
        Ar31 = Ar3;
        Ar32 = Ar3 + Na;
    else
        Ar32 = Ar3;
        Ar31 = Ar3 + Na;
    end
    Tri_ar2(4*K-3,:) = [Ar11 2*Na+K*3 Ar32];
    Tri_ar2(4*K-2,:) = [Ar12 Ar21 2*Na+K*3-2];
    Tri_ar2(4*K-1,:) = [Ar22 Ar31 2*Na+K*3-1];
    Tri_ar2(4*K,:) = [2*Na+K*3-2 2*Na+K*3-1 2*Na+K*3];
    if Ar_tri2(Ar11,1) == 0
        Ar_tri2(Ar11,1) = 4*K-3;
    else
        Ar_tri2(Ar11,2) = 4*K-3;
    end
    if Ar_tri2(Ar12,1) == 0
        Ar_tri2(Ar12,1) = 4*K-2;
    else
        Ar_tri2(Ar12,2) = 4*K-2;
    end
    if Ar_tri2(Ar31,1) == 0
        Ar_tri2(Ar31,1) = 4*K-1;
    else
        Ar_tri2(Ar31,2) = 4*K-1;
    end
    if Ar_tri2(Ar32,1) == 0
        Ar_tri2(Ar32,1) = 4*K-3;
    else
        Ar_tri2(Ar32,2) = 4*K-3;
    end
    if Ar_tri2(Ar22,1) == 0
        Ar_tri2(Ar22,1) = 4*K-1;
    else
        Ar_tri2(Ar22,2) = 4*K-1;
    end
    if Ar_tri2(Ar21,1) == 0
        Ar_tri2(Ar21,1) = 4*K-2;
    else
        Ar_tri2(Ar21,2) = 4*K-2;
    end
    Ar_tri2(2*Na+K*3-2,:) = [4*K-2 4*K];
    Ar_tri2(2*Na+K*3-1,:) = [4*K-1 4*K];
    Ar_tri2(2*Na+K*3,:) = [4*K-3 4*K];
end
Barycentres2 = zeros(Na2,2);
for F=1:Na2
    Barycentres2(F,:) = (Coorneu2(Numaretes2(F,1),:) + Coorneu2(Numaretes2(F,2),:))/2;
end
%{
figure;

% maillage
trimesh(Numtri2(:,1:3),Coorneu2(:,1),Coorneu2(:,2),zeros(Nbpt2,1));
hold on;
plot(Coorneu2(:,1),Coorneu2(:,2),'.r', "markersize", 10);
view(2);
%}