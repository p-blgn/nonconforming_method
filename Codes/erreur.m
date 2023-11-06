nom_maillage = "geom_non_homogene.msh";
nombre_de_raffinements = 4;
S = zeros(nombre_de_raffinements+1,1);
H = zeros(nombre_de_raffinements+1,1);
E = zeros(nombre_de_raffinements+1,1);
Dim = zeros(nombre_de_raffinements+1,1);
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes_ext,Numaretes_bords] = lecture_msh(nom_maillage);
[Barycentres,Na,Numaretes,isavertex,Ar_tri,Nbaretes_int,Tri_ar] = Initialisation_matrices(Coorneu,Nbtri,Numtri,Nbaretes_ext,Numaretes_bords,Nbpt);
[U_h,Nbtri,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,KK] = test_non_h(Barycentres,Tri_ar,Ar_tri,Nbaretes_int,Na,Numaretes,Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Numaretes_bords);
%[U_h,Nbtri,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,KK] = test(Ar_tri,Nbaretes_int,Na,Numaretes,Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Numaretes_bords);
Aire_tri = Aire_triangles(Coorneu,Nbtri,Numtri);
longueurs_faces = longueurs(Na,Coorneu,Numaretes);
G_U = zeros(Na,2); %solution exacte
%G_U(:,1) = pi*cos(pi*Barycentres(:,1)).*sin(pi*Barycentres(:,2));
%G_U(:,2) = pi*sin(pi*Barycentres(:,1)).*cos(pi*Barycentres(:,2));

for i=1:Na
    G_U(i,:) = grad_polaire(Barycentres(i,1),Barycentres(i,2))';
end

%Val_sommets = vers_sommets_non_homogene(U_h,Nbtri,Tri_ar,Na,Numaretes);
%Val_sommets2 = reconstruction_potentiel_non_homogene(Val_sommets,Numtri,Ns,Nbtri,Refneu,Coorneu);
%Val_sommets = vers_sommets(U_h,Nbtri,Tri_ar,Na,Numaretes);
%Val_sommets2 = reconstruction_potentiel(Val_sommets,Numtri,Ns,Nbtri,Refneu);
%{
U_h = zeros(Na,1); %valeur de s_h sur les barycentres, on peut alors appliquer la formule pour avoir le gradient
for F=1:Na
    U_h(F) = (Val_sommets2(Numaretes(F,1)) + Val_sommets2(Numaretes(F,2)))/2;
end
%}
s = 0;
for K=1:Nbtri
    a = 0;
    for F=Tri_ar(K,:) %Calculer vec une fois plutôt que 3
        vec = zeros(2,1);
        %Ici on utilise le fait que le vecteur normal est le vecteur
        %[arête1 arête2] après une rotation de -pi/2 car on tourne en sens
        %trigo, si on tourne dans le sens horaire on aura un signe d'écart
        vec(1) = (U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),2)-Coorneu(Numtri(K,1),2)) + U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),2)-Coorneu(Numtri(K,2),2)) + U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),2)-Coorneu(Numtri(K,3),2)))/Aire_tri(K);
        vec(2) = (-U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),1)-Coorneu(Numtri(K,1),1)) - U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),1)-Coorneu(Numtri(K,2),1)) - U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),1)-Coorneu(Numtri(K,3),1)))/Aire_tri(K);
        a = a + (G_U(F,:)-vec')*(G_U(F,:)'-vec);
    end
    s = s + a*Aire_tri(K)/3;
end
E(1) = sqrt(sum(estimation(U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,Barycentres,Aire_tri)));
S(1) = sqrt(s);
Dim(1) = Na;
H(1) = max(longueurs_faces);
for i=1:nombre_de_raffinements
    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Numaretes_bords,Na,Barycentres,Numaretes,Ar_tri,Nbaretes_int,Tri_ar] = raffinement(Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Na,Barycentres,Numaretes,Nbaretes_int,Tri_ar);
    [U_h,Nbtri,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,KK] = test_non_h(Barycentres,Tri_ar,Ar_tri,Nbaretes_int,Na,Numaretes,Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Numaretes_bords);
    %[U_h,Nbtri,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,KK] = test(Ar_tri,Nbaretes_int,Na,Numaretes,Nbpt,Nbtri,Coorneu,Refneu,Numtri,Nbaretes_ext,Numaretes_bords);
    Aire_tri = Aire_triangles(Coorneu,Nbtri,Numtri);
    longueurs_faces = longueurs(Na,Coorneu,Numaretes);
    G_U = zeros(Na,2);
    %G_U(:,1) = pi*cos(pi*Barycentres(:,1)).*sin(pi*Barycentres(:,2));
    %G_U(:,2) = pi*sin(pi*Barycentres(:,1)).*cos(pi*Barycentres(:,2));
    
    for m=1:Na
        G_U(m,:) = grad_polaire(Barycentres(m,1),Barycentres(m,2))';
    end
    
    %Val_sommets = vers_sommets_non_homogene(U_h,Nbtri,Tri_ar,Na,Numaretes);
    %Val_sommets2 = reconstruction_potentiel_non_homogene(Val_sommets,Numtri,Ns,Nbtri,Refneu,Coorneu);
    %Val_sommets = vers_sommets(U_h,Nbtri,Tri_ar,Na,Numaretes);
    %Val_sommets2 = reconstruction_potentiel(Val_sommets,Numtri,Ns,Nbtri,Refneu);
    %{
    U_h = zeros(Na,1); %valeur de s_h sur les barycentres, on peut alors appliquer la formule pour avoir le gradient
    for F=1:Na
        U_h(F) = (Val_sommets2(Numaretes(F,1)) + Val_sommets2(Numaretes(F,2)))/2;
    end
    %}
    s = 0;
    erreurs = zeros(Nbtri,1);
    for K=1:Nbtri
        a = 0;
        for F=Tri_ar(K,:)
            vec = zeros(2,1);
            vec(1) = (U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),2)-Coorneu(Numtri(K,1),2)) + U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),2)-Coorneu(Numtri(K,2),2)) + U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),2)-Coorneu(Numtri(K,3),2)))/Aire_tri(K);
            vec(2) = (-U_h(Tri_ar(K,1))*(Coorneu(Numtri(K,2),1)-Coorneu(Numtri(K,1),1)) - U_h(Tri_ar(K,2))*(Coorneu(Numtri(K,3),1)-Coorneu(Numtri(K,2),1)) - U_h(Tri_ar(K,3))*(Coorneu(Numtri(K,1),1)-Coorneu(Numtri(K,3),1)))/Aire_tri(K);
            a = a + (G_U(F,:)-vec')*(G_U(F,:)'-vec);
        end
        erreurs(K) = a*Aire_tri(K)/3;
        s = s + a*Aire_tri(K)/3;
    end
    estimations = estimation(U_h,Nbtri,Tri_ar,Na,Numaretes,Numtri,Coorneu,Ns,Refneu,Barycentres,Aire_tri);
    E(i+1) = sqrt(sum(estimations));
    S(i+1) = sqrt(s);
    Dim(i+1) = Na;
    H(i+1) = max(longueurs_faces);
end


%figure;
%{
p = polyfit(log(H),log(S),1);
y1 = H.^(p(1))*exp(p(2));
loglog(H,S,'b-*');
hold on;
loglog(H,y1,'r--');
hold on;
legend("",'Location','southwest')
%}


%p = polyfit(log(Dim),log(S),1);
%y1 = Dim.^(p(1))*exp(p(2));
%loglog(Dim,S,'b-*');
%hold on;
%loglog(Dim,E,'g-*');
%hold on;
%legend("erreur","estimateur");
%xlabel('Nombre de faces');
%loglog(Dim,y1,'r--');
%hold on;

figure;
semilogx(Dim,E./S);
xlabel('Nombre de faces');
ylabel("Indice d'efficacité");


figure;
for K=1:Nbtri
    S1 = Numtri(K,1);
    S2 = Numtri(K,2);
    S3 = Numtri(K,3);
    T = [1 2 3];
    X = [Coorneu(S1,1);Coorneu(S2,1);Coorneu(S3,1)];
    Y = [Coorneu(S1,2);Coorneu(S2,2);Coorneu(S3,2)];
    Z = [sqrt(estimations(K));sqrt(estimations(K));sqrt(estimations(K))];
    trisurf(T,X,Y,Z);
    hold on;
end
view(2);
colorbar;
figure;
for K=1:Nbtri
    S1 = Numtri(K,1);
    S2 = Numtri(K,2);
    S3 = Numtri(K,3);
    T = [1 2 3];
    X = [Coorneu(S1,1);Coorneu(S2,1);Coorneu(S3,1)];
    Y = [Coorneu(S1,2);Coorneu(S2,2);Coorneu(S3,2)];
    Z = [sqrt(erreurs(K));sqrt(erreurs(K));sqrt(erreurs(K))];
    trisurf(T,X,Y,Z);
    hold on;
end
view(2);
colorbar;
