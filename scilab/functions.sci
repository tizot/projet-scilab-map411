///////////////////////////////
// Un calcul de flot optique //
//      BRUNEAU Basile       //
//      MASSET Camille       //
//          MAP411           //
///////////////////////////////

// Bibliothèque de fonctions utilisées dans 'projet.sce'

// QUESTION 4
// Approximation des dérivées en x, y et t de l'intensité du flot optique
function [derX, derY, derT] = derivees(I1, I2)
    if (~issquare(I1) | ~issquare(I2)) then
        error("L''une des matrices n''est pas carrée.")
    end
    if (length(I1(1,:)) <> length(I2(1,:))) then
        error("Les matrices ne sont pas de même dimension.")
    end
    
    N = length(I1(1,:)); // dimension des matrices
    derX = zeros(I1); // initialisation du résultat
    derY = zeros(I1);
    derT = zeros(I1);
    
    // dérivées intérieures
    derX(1:N-1, 1:N-1) = (I1(1:N-1, 2:N) - I1(1:N-1, 1:N-1) + I1(2:N, 2:N) - I1(2:N, 1:N-1) + I2(1:N-1, 2:N) - I2(1:N-1, 1:N-1) + I2(2:N, 2:N) - I2(2:N, 1:N-1))/4;
    derY(1:N-1, 1:N-1) = (I1(2:N, 1:N-1) - I1(1:N-1, 1:N-1) + I1(2:N, 2:N) - I1(1:N-1, 2:N) + I2(2:N, 1:N-1) - I2(1:N-1, 1:N-1) + I2(2:N, 2:N) - I2(1:N-1, 2:N))/4;
    derT(1:N-1, 1:N-1) = (I2(1:N-1, 1:N-1) - I1(1:N-1, 1:N-1) + I2(2:N, 1:N-1) - I1(2:N, 1:N-1) + I2(2:N, 2:N) - I1(2:N, 2:N) + I2(1:N-1, 2:N) - I1(1:N-1, 2:N))/4;
    
    // dérivées au bord (condition de Neumann)
    // à vérifier
    derX(N, 1:N-1) = (I1(N, 2:N) - I1(N, 1:N-1) + 0 + I2(N, 2:N) - I2(N, 1:N-1) + 0)/4;
    // derX(1:N-1, N) = (0 + 0 + 0 + 0)/4;
    // derX(N, N) = 0;
    
    // derY(N, 1:N-1) = (0 + 0 + 0 + 0)/4;
    derY(1:N-1, N) = (I1(2:N, N) - I1(1:N-1, N) + 0 + I2(2:N, N) - I2(1:N-1, N) + 0)/4;
    // derY(N, N) = 0;
    
    derT(N, 1:N-1) = (I2(N, 1:N-1) - I1(N, 1:N-1) + I2(N, 2:N) - I1(N, 2:N))/2;
    derT(1:N-1, N) = (I2(1:N-1, N) - I1(1:N-1, N) + I2(2:N, N) - I1(2:N, N))/2;
    derT(N, N) = I2(N, N) - I1(N, N);
endfunction

// QUESTION 5
// Calcul d'une matrice moyenne
function Umoy = moyenneMat(U)
    [n,p] = size(U); // U a n lignes et p colonnes (on supposera n,p > 2)
    Umoy = zeros(U);
    
    // on pondère les éléments qui ont un côté en commun 2 fois plus que les éléments qui n'ont qu'un coin en commun
    // éléments intérieurs
    Umoy(2:n-1, 2:p-1) = 1/6 * (U(1:n-2, 2:p-1) + U(3:n, 2:p-1) + U(2:n-1, 3:p) + U(2:n-1, 1:p-2)) + 1/12 * (U(1:n-2, 1:p-2) + U(1:n-2, 3:p) + U(3:n, 1:p-2) + U(3:n, 3:p));
    
    // éléments aux bords
    Umoy(1, 2:p-1) = 1/4 * (U(1, 1:p-2) + U(1, 3:p) + U(2, 2:p-1)) + 1/8 * (U(2, 1:p-2) + U(2, 3:p)); // ligne 1
    Umoy(n, 2:p-1) = 1/4 * (U(n, 1:p-2) + U(n, 3:p) + U(n-1, 2:p-1)) + 1/8 * (U(n-1, 1:p-2) + U(n-1, 3:p)); // ligne n
    Umoy(2:n-1, 1) = 1/4 * (U(1:n-2, 1) + U(3:n, 1) + U(2:n-1, 2)) + 1/8 * (U(1:n-2, 2) + U(3:n, 2)); // colonne 1
    Umoy(2:n-1, p) = 1/4 * (U(1:n-2, p) + U(3:n, p) + U(2:n-1, p-1)) + 1/8 * (U(1:n-2, p-1) + U(3:n, p-1)); // colonne p
    
    // éléments aux coins
    Umoy(1, 1) = 2/5 * (U(1, 2) + U(2, 1)) + 1/5 * U(2, 2);
    Umoy(1, p) = 2/5 * (U(1, p-1) + U(2, p)) + 1/5 * U(2, p-1);
    Umoy(n, 1) = 2/5 * (U(n-1, 1) + U(n, 2)) + 1/5 * U(n-1, 2);
    Umoy(n, p) = 2/5 * (U(n-1, p) + U(n, p-1)) + 1/5 * U(n-1, p-1);
endfunction
