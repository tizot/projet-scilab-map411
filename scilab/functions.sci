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
    derX(1:N-1, N) = (0 + 0 + 0 + 0)/4;
    derX(N, N) = 0;
    
    derY(N, 1:N-1) = (0 + 0 + 0 + 0)/4;
    derY(1:N-1, N) = (I1(2:N, N) - I1(1:N-1, N) + 0 + I2(2:N, N) - I2(1:N-1, N) + 0)/4;
    derY(N, N) = 0;
    
    derT(N, 1:N-1) = (I2(N, 1:N-1) - I1(N, 1:N-1) + 0 + 0 + I2(N, 2:N) - I1(N, 2:N))/4;
    derT(1:N-1, N) = (I2(1:N-1, N) - I1(1:N-1, N) + I2(2:N, N) - I1(2:N, N) + 0 + 0)/4;
    derT(N, N) = 0;
endfunction
