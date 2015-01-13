///////////////////////////////
// Un calcul de flot optique //
//      BRUNEAU Basile       //
//      MASSET Camille       //
//          MAP411           //
///////////////////////////////

// Programme Scilab du projet.
// Les fonctions utilisées ont été définies dans 'functions.sci'

exec('functions.sci', -1);

// ===== On applique Q7
// Paramètres
niter = 64; // nombre maximal d'itérations
epsilon = 0.001; // erreur tolérée
alpha = 0.1;

I1 = fscanfMat("I1.txt");
I2 = fscanfMat("I2.txt");
[N, M] = size(I1)

[U, V] = flow(I1, I2, niter, epsilon, alpha)
clf()
champ(1:N/5, 1:M/5,  V(1:5:N, 1:5:M), U(1:5:N,1:5:M)) // , arfact=0.1)
sleep(4000)


// Et maintenant on applique Q14
[U, V] = smartFlow(I1, I2)
clf()
champ(1:N/5, 1:M/5, V(1:5:N, 1:5:M), U(1:5:N,1:5:M))

