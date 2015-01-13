///////////////////////////////
// Un calcul de flot optique //
//      BRUNEAU Basile       //
//      MASSET Camille       //
//          MAP411           //
///////////////////////////////

// Programme Scilab du projet.
// Les fonctions utilisées ont été définies dans 'functions.sci'

exec('functions.sci', -1);

// clf()
// I1 = fscanfMat("I1.txt");
// G = blur(I1, 8, 0.05)
// toGray(G)
// fprintfMat('out/grove-blur-8-0.05.txt', G)
// return

// sleep(2000)

// Question 7
// Implémentation de la méthode de Horn et Schunk

// Paramètres
niter = 64; // nombre maximal d'itérations
epsilon = 0.001; // erreur tolérée
alpha = 1;

I1 = fscanfMat("I1.txt");
I2 = fscanfMat("I2.txt");
[N, M] = size(I1)

[U, V] = flow(I1, I2, niter, epsilon, alpha)
clf()
champ(1:N/5, 1:M/5, U(1:5:N,1:5:M), V(1:5:N, 1:5:M)) // , arfact=0.1)
sleep(4000)

// I1 = blur(I1, 2, 0.05)
// I2 = blur(I2, 2, 0.05)

// [U, V] = flow(I1, I2, niter, epsilon, alpha)

[U, V] = smartFlow(I1, I2)
clf()
champ(1:N/5, 1:M/5, U(1:5:N,1:5:M), V(1:5:N, 1:5:M))

//G = convolKernel(2, 0.05)

//toGray(blur(I1, 20, 0.05))
