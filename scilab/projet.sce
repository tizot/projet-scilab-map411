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
niter = 100; // nombre maximal d'itérations
epsilon = 0.01; // erreur tolérée
alpha = 1;

I1 = RGB2Gray(ReadImage('frame10.png'));
I2 = RGB2Gray(ReadImage('frame11.png'));

// Valeurs initiales
err = 2*epsilon;
i = 1;
U = zeros(I1);
V = zeros(I1);

// Dérivées de I
[Ix, Iy, It] = derivees(I1, I2);

while (i <= niter & err > epsilon)
    i++;
    [Up, Vp] = iter(U, V, Ix, Iy, It, alpha);
    err = trace((Up-U)' * (Up-U)) + trace((Vp-V)' * (Vp-V));
    U = Up;
    V = Vp;
end

