function Longueurs_Ar = longueurs(Na,Coorneu,Numaretes)
Longueurs_Ar = zeros(Na,1);
for i=1:Na
    sommet1 = Numaretes(i,1);
    sommet2 = Numaretes(i,2);
    x1 = Coorneu(sommet1,1);
    x2 = Coorneu(sommet2,1);
    y1 = Coorneu(sommet1,2);
    y2 = Coorneu(sommet2,2);
    Longueurs_Ar(i) = norm([x2 - x1;y2 - y1]);
end