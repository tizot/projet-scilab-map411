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
    [n,p] = size(I1); // I1 a n lignes et p colonnes (on supposera n,p > 2)
    
    if (size(I2) <> [n,p]) then
        error("Les matrices ne sont pas de même dimension.")
    end
    
    derX = zeros(I1); // initialisation du résultat
    derY = zeros(I1);
    derT = zeros(I1);
    
    // dérivées intérieures
    derX(1:n-1, 1:p-1) = (I1(1:n-1, 2:p) - I1(1:n-1, 1:p-1) + I1(2:n, 2:p) - I1(2:n, 1:p-1) + I2(1:n-1, 2:p) - I2(1:n-1, 1:p-1) + I2(2:n, 2:p) - I2(2:n, 1:p-1))/4;
    derY(1:n-1, 1:p-1) = (I1(2:n, 1:p-1) - I1(1:n-1, 1:p-1) + I1(2:n, 2:p) - I1(1:n-1, 2:p) + I2(2:n, 1:p-1) - I2(1:n-1, 1:p-1) + I2(2:n, 2:p) - I2(1:n-1, 2:p))/4;
    derT(1:n-1, 1:p-1) = (I2(1:n-1, 1:p-1) - I1(1:n-1, 1:p-1) + I2(2:n, 1:p-1) - I1(2:n, 1:p-1) + I2(2:n, 2:p) - I1(2:n, 2:p) + I2(1:n-1, 2:p) - I1(1:n-1, 2:p))/4;
    
    // dérivées au bord (condition de Neumann)
    // à vérifier
    derX(n, 1:p-1) = (I1(n, 2:p) - I1(n, 1:p-1) + 0 + I2(n, 2:p) - I2(n, 1:p-1) + 0)/4;
    // derX(1:N-1, N) = (0 + 0 + 0 + 0)/4;
    // derX(N, N) = 0;
    
    // derY(N, 1:N-1) = (0 + 0 + 0 + 0)/4;
    derY(1:n-1, p) = (I1(2:n, p) - I1(1:n-1, p) + 0 + I2(2:n, p) - I2(1:n-1, p) + 0)/4;
    // derY(N, N) = 0;
    
    derT(n, 1:p-1) = (I2(n, 1:p-1) - I1(n, 1:p-1) + I2(n, 2:p) - I1(n, 2:p))/4;
    derT(1:n-1, p) = (I2(1:n-1, p) - I1(1:n-1, p) + I2(2:n, p) - I1(2:n, p))/4;
    derT(n, p) = (I2(n, p) - I1(n, p))/4;
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
    Umoy(1, 2:p-1) = 1/6 * (U(1, 1:p-2) + U(1, 3:p) + U(2, 2:p-1) + U(1, 2:p-1)) + 1/12 * (U(2, 1:p-2) + U(2, 3:p) + U(1, 1:p-2) + U(1, 3:p)); // ligne 1
    Umoy(n, 2:p-1) = 1/6 * (U(n, 1:p-2) + U(n, 3:p) + U(n-1, 2:p-1) + U(n, 2:p-1)) + 1/12 * (U(n-1, 1:p-2) + U(n-1, 3:p) + U(n, 1:p-2) + U(n, 3:p)); // ligne n
    Umoy(2:n-1, 1) = 1/6 * (U(1:n-2, 1) + U(3:n, 1) + U(2:n-1, 2) + U(2:n-1, 1)) + 1/12 * (U(1:n-2, 2) + U(3:n, 2) + U(1:n-2, 1) + U(3:n, 1)); // colonne 1
    Umoy(2:n-1, p) = 1/6 * (U(1:n-2, p) + U(3:n, p) + U(2:n-1, p-1) + U(2:n-1, p)) + 1/12 * (U(1:n-2, p-1) + U(3:n, p-1) + U(1:n-2, p) + U(3:n, p)); // colonne p
    
    // éléments aux coins
    // A CHANGER - TODO
    Umoy(1, 1) = 2/5 * (U(1, 2) + U(2, 1)) + 1/5 * U(2, 2);
    Umoy(1, p) = 2/5 * (U(1, p-1) + U(2, p)) + 1/5 * U(2, p-1);
    Umoy(n, 1) = 2/5 * (U(n-1, 1) + U(n, 2)) + 1/5 * U(n-1, 2);
    Umoy(n, p) = 2/5 * (U(n-1, p) + U(n, p-1)) + 1/5 * U(n-1, p-1);
endfunction

function [Unp, Vnp] = iter(Un, Vn, Ix, Iy, It, alpha)
    UnMoy = moyenneMat(Un);
    VnMoy = moyenneMat(Vn);
    matAlpha2 = alpha^2 * ones(Un);
    frac = (Ix .* UnMoy + Iy .* VnMoy + It) .* ((matAlpha2 + Ix .* Ix + Iy .* Iy).^(-1));
    Unp = UnMoy - Ix .* frac;
    Vnp = VnMoy - Iy .* frac;
endfunction

function G = convolKernel(sigma, eta)
    // Indice limite
    ilim = ceil(sqrt(-2 * sigma^2 * log(eta)));
    
    G = zeros(2 * ilim + 1);
    
    // Remplissage de la matrice
    for i=-ilim:ilim
        for j=-ilim:ilim
            G(i+ilim+1,j+ilim+1) = exp(-(i^2 + j^2)/(2 * sigma^2))/(2 * %pi * sigma^2);
        end
    end
endfunction

function refI = reflectMat(I, n)
    [N, M] = size(I); // on suppose n < N,M
    refI = zeros(N + 2*n, M + 2*n);
    refI(n+1:N+n, n+1:M+n) = I;
   
    for i=1:n
        refI(n+1-i, n+1:M+n) = I(i, :); // lignes "au-dessus"
        refI(N+n+1+i, n+1:M+n) = I(N-i, :); // lignes "en-dessous"
        refI(n+1:N+n, n+1-i) = I(:, i); // colonnes "à gauche"
        refI(n+1:N+n, M+n+1+i) = I(:, M-i); // colonnes "à droite"
    end
endfunction

function B = blur(I, sigma, eta)
    G = convolKernel(sigma, eta);
    
    n = length(G(1, :));
    refI = reflectMat(I, n);
    
    C = conv2(refI, G, "same");
    
    [N, M] = size(I);
    B = C(n+1:N+n, n+1:M+n);
endfunction

function toGray(I)
    [N, M] = size(I);
    cmap = graycolormap(256);
    grayplot(1:N, 1:M, I)
    xset('colormap', cmap);
endfunction
