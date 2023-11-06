function Aire_tri = Aire_triangles(Coorneu,Nbtri,Numtri)
Aire_tri = zeros(Nbtri,1);
for l=1:Nbtri
    sommet1 = Coorneu(Numtri(l,1),:);
    sommet2 = Coorneu(Numtri(l,2),:);
    sommet3 = Coorneu(Numtri(l,3),:);
    x1 = sommet2(1) - sommet1(1);
    x2 = sommet3(1) - sommet1(1);
    y1 = sommet2(2) - sommet1(2);
    y2 = sommet3(2) - sommet1(2);
    Aire_tri(l) = (x1*y2 - y1*x2)/2;
end