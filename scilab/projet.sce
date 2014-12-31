///////////////////////////////
// Un calcul de flot optique //
//      BRUNEAU Basile       //
//      MASSET Camille       //
//          MAP411           //
///////////////////////////////

// Programme Scilab du projet.
// Les fonctions utilisées ont été définies dans 'functions.sci'

exec('functions.sci', -1);

I1 = rand(10, 10, "normal");
I2 = rand(10, 10, "normal");

Umoy = moyenneMat(I1);
I1
