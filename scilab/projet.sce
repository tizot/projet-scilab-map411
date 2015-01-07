///////////////////////////////
// Un calcul de flot optique //
//      BRUNEAU Basile       //
//      MASSET Camille       //
//          MAP411           //
///////////////////////////////

// Programme Scilab du projet.
// Les fonctions utilisées ont été définies dans 'functions.sci'

exec('functions.sci', -1);


// Question 7
// Implémentation de la méthode de Horn et Schunk

// Paramètres
niter = 640; // nombre maximal d'itérations
epsilon = 0.001; // erreur tolérée
alpha = 0.01;

I1 = fscanfMat("I1.txt");
I2 = fscanfMat("I2.txt");

[N, M] = size(I1)

// Valeurs initiales
err = 1;
i = 1;
U = zeros(I1);
V = zeros(I1);

// Dérivées de I
[Ix, Iy, It] = derivees(I1, I2);

while (i <= niter & err > epsilon)
    i = i + 1
    [Up, Vp] = iter(U, V, Ix, Iy, It, alpha);
    err = trace((Up-U)' * (Up-U)) + trace((Vp-V)' * (Vp-V))
    U = Up;
    V = Vp;
    //champ(1:N, 1:M, round(U*10)/10, round(V*10)/10)
    //sleep(2000)
    //clf()
end

disp(i)
disp(err)

clf()
champ(1:N/5, 1:M/5, U(1:5:N,1:5:M), V(1:5:N, 1:5:M))

G = convolKernel(2, 0.05)
