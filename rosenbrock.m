% x : valeur de la variable d'optimisation
% params : structure contenant les parametres necessaires pour evaluer la fonction objectif
% f,g,h : valeur de la fonction objectif, son gradient et son hessien
%

% A COMPLETER
function [f, g, H, J] = rosenbrock(x, params)
    n = numel(x);
    %calcul de f
    f = 0;
    for i=1:n-1
        f = f + params*(x(i+1) - x(i)^2)^2 + (1 - x(i))^2;
    end
    %calcul de g
    g = zeros(n,1);
    g(1) = -4*params*x(1)*(x(2)-x(1)^2)-2*(1-x(1));
    g(n) = 2*params*(x(n)-x(n-1)^2);
    if n > 2 
        for i =2:n-1
            g(i) = -4*params*x(i)*(x(i+1)-x(i)^2)-2*(1-x(i)) + 2*params*(x(i)-x(i-1)^2);
        end
    end
    %calcul de H
    H = zeros(n,n);
    %calcul des éléments diagonaux
    H(1,1) = 2-4*params*(x(2)-3*x(1)^2);
    H(n,n) = 2*params;
    for i=2:n-1
        H(i,i) = 2-4*params*(x(i+1)-3*x(i)^2) + 2*params;
    end
    %calcul des éléments de la diagonale supérieure et inférieure 
    for i=1:n
        for j=2:n
            if i == j-1
                H(i,j) = -4*params*x(i);
                H(j,i) = -4*params*x(i);
            end
        end
    end
    %calcul de la jacobienne
    J = zeros(2*(n-1), n);
    for i = 1:n-1
        %première ligne du jacobien pour r1(x_i,x_{i+1})
        J(2*i, i) =  -2*sqrt(params)*x(i);
        J(2*i, i+1) = sqrt(params);
        %deuxième ligne du jacobien pour r(x_i)
        J(2*i-1,i) = -1;
    end
          
end



    
