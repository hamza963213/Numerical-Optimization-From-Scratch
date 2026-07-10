% critfun : nom du fichier .m ?valuant la fonction objectif, son gradient et hessien
% params : param?tres pour l'?valuation de la fonction objectif
% options : options n?cessaires pour la mise en oeuvre des algorithmes
%    options.method : 'gradient', 'gradient conjuge', 'Newton', 'Quasi-Newton', ...
%    options.pas : 'fixe', 'variable'
%    options.const : constante d'armijo
%    options.beta : taux de rebroussement
%    options.tolX, options.tolF, options.tolG, options.maxiter
% x0 : point initial
% xh : point final
% result : structure contenant les r?sultats
%    result.iter : nombre d'iterations r?alis?es
%    result.crit : valeur du crit?re par iteration
%    result.grad : norm du gradient par iteration 
%    result.temp : temps de calcul
%    result.stop : condition d'arret ('TolX', 'TolG', 'TolF', 'Maxiter')
% xval : valeurs des itérées
% A COMPLETER

function [xh,result,xval] = optimdescent(critfun,params,options,x0)
    %démarrage du chronomètre pour mesurer le temps de calcul
    tic;
    %acquistion de la valeur du pas
    if strcmp(options.pas,'fixe')
        pas= input('donner la valeur du pas fixe :');
    else 
        pas = input('donner la valeur du pas initial :');
    end
    %initialisation des itérées/calcul de la première itération
    [f_prec,g_prec,H_prec] = critfun(x0,params);
    x_prec= x0;
    d = -g_prec;
    xh = x_prec + pas * d;
    [f_h,g_h,H_h, J_h] = rosenbrock(xh,params);
    
    %initialisation de B_h pour Quasi-newton
    B_h = H_h;
    %choix de lambda si Levenberg-Marquardt
    if strcmp(options.method, 'Levenberg-Marquardt')
        lambda = input('choisir lambda :');
    end
    %initialisation des variables résultats
    xval = [x0, xh];
    result.crit = [f_prec, f_h];
    result.grad= [norm(g_prec) , norm(g_h)];
    k = 1;
    while k < options.maxiter && norm(g_h) > options.tolG && norm((xh-x_prec)/xh) > options.tolX && norm((f_h-f_prec)/f_h) > options.tolF
        switch options.method
            case 'gradient'
                %choix de la direction d
                d = -g_h;
                %vérifier que d est bien une descente
                if g_h'*d > 0 
                    d=-d;
                end
                %sauvegarde itération précédente pour les conditions d'arrêt
                g_prec = g_h;
                x_prec = xh;
                f_prec = f_h;
                %calcul de l'itération en cours 
                xh = x_prec + pas * d;
                [f_h,g_h,~] = rosenbrock(xh,params);
                %backtracking
                if strcmp(options.pas,'variable')
                    while (f_h-f_prec) > 1e-4*pas*g_prec*d'
                        pas = 0.75*pas;
                        xh = x_prec + pas * d;
                        [f_h,g_h,~] = rosenbrock(xh,params);
                    end
                end
                %sauvegarde des résultats intermédiares
                result.crit = [result.crit, f_h]; %critère par itération 
                result.grad = [result.grad, norm(g_h)]; %gradient par itération
                xval = [xval, xh]; %sauvegarde des itérées
            case 'gradient conjuge'       
                %calcul de la direction d
                betha_h = norm(g_h)^2/norm(g_prec)^2;
                d = -g_h + betha_h * d;
                %vérifier que d est bien une descente
                if g_h'*d > 0 
                    d=-d;
                end
                %sauvegarde itération précédente pour les conditions d'arrêt
                g_prec = g_h;
                x_prec = xh;
                f_prec = f_h;
                %calcul de l'itération en cours/recherche d'un bon pas
                xh = x_prec + pas * d;
                [f_h,g_h,~] = rosenbrock(xh,params);
                %backtracking
                if strcmp(options.pas,'variable')
                    while (f_h-f_prec) > 1e-4*pas*g_prec*d'
                        pas = 0.75*pas;
                        xh = x_prec + pas * d;
                        [f_h,g_h,~] = rosenbrock(xh,params);
                    end
                end
                %sauvegarde des résultats intermédiares
                result.crit = [result.crit, f_h]; %critère par itération 
                result.grad = [result.grad, norm(g_h)]; %gradient par itération
                xval = [xval, xh]; %sauvegarde des itérées

            case 'Newton'
                %calcul de la direction d
                d = linsolve(H_h,-g_h);
         
                %vérifier que d est bien une descente
                if g_h'*d > 0 
                    d=-d;
                end
                %sauvegarde itération précédente pour les conditions d'arrêt
                g_prec = g_h;
                x_prec = xh;
                f_prec = f_h;
                %calcul de l'itération en cours/recherhce d'un bon pas
                xh = x_prec + pas * d;
                [f_h,g_h,H_h] = rosenbrock(xh,params);
                %backtracking
                if strcmp(options.pas,'variable')
                    while (f_h-f_prec) > 1e-4*pas*g_prec*d'
                        pas = 0.75*pas;
                        xh = x_prec + pas * d;
                        [f_h,g_h,H_h] = rosenbrock(xh,params);
                    end
                end
                %sauvegarde des résultats intermédiares
                result.crit = [result.crit, f_h]; %critère par itération 
                result.grad = [result.grad, norm(g_h)]; %gradient par itération
                xval = [xval, xh]; %sauvegarde des itérées

            case 'Quasi-Newton'
                s = xh - x_prec;
                y = g_h -g_prec;
                %calcul de la direction d    
                B_h = B_h - (B_h*s*(s)'*B_h)/(s'*B_h*s) + (y*y')/(y'*s);
                d = linsolve(B_h,-g_h);
                
                %vérifier que d est bien une descente
                if g_h'*d > 0 
                    d=-d;
                end
                %sauvegarde itération précédente pour les conditions d'arrêt
                g_prec = g_h;
                x_prec = xh;
                f_prec = f_h;
                %calcul de l'itération en cours/recherhce d'un bon pas
                xh = x_prec + pas * d;
                [f_h,g_h,~] = rosenbrock(xh,params);
                %backtracking
                if strcmp(options.pas,'variable')
                    while (f_h-f_prec) > 1e-4*pas*g_prec*d'
                        pas = 0.75*pas;
                        xh = x_prec + pas * d;
                        [f_h,g_h,~] = rosenbrock(xh,params);
                    end
                end
                %sauvegarde des résultats intermédiares
                result.crit = [result.crit, f_h]; %critère par itération 
                result.grad = [result.grad, norm(g_h)]; %gradient par itération
                xval = [xval, xh]; %sauvegarde des itérées
            case 'Gauss-Newton'
                %calcul de la direction d
                d = linsolve((J_h)'*J_h, -g_h);
                 %vérifier que d est bien une descente
                if g_h'*d > 0 
                    d=-d;
                end
                %sauvegarde itération précédente pour les conditions d'arrêt
                g_prec = g_h;
                x_prec = xh;
                f_prec = f_h;
                %calcul de l'itération en cours/recherhce d'un bon pas
                xh = x_prec + pas * d;
                [f_h,g_h,~, J_h] = rosenbrock(xh,params);
                %backtracking
                if strcmp(options.pas,'variable')
                    while (f_h-f_prec) > 1e-4*pas*g_prec*d'
                        pas = 0.75*pas;
                        xh = x_prec + pas * d;
                        [f_h,g_h,~, J_h] = rosenbrock(xh,params);
                    end
                end
                %sauvegarde des résultats intermédiares
                result.crit = [result.crit, f_h]; %critère par itération 
                result.grad = [result.grad, norm(g_h)]; %gradient par itération
                xval = [xval, xh]; %sauvegarde des itérées
            case 'Levenberg-Marquardt'
                d = linsolve((J_h)'*J_h + lambda*diag((J_h)'*J_h), -g_h);
                 %vérifier que d est bien une descente
                if g_h'*d > 0 
                    d=-d;
                end
                %sauvegarde itération précédente pour les conditions d'arrêt
                g_prec = g_h;
                x_prec = xh;
                f_prec = f_h;
                %calcul de l'itération en cours/recherhce d'un bon pas
                xh = x_prec + pas * d;
                [f_h,g_h,~, J_h] = rosenbrock(xh,params);
                %backtracking
                if strcmp(options.pas,'variable')
                    while (f_h-f_prec) > 1e-4*pas*g_prec*d'
                        pas = 0.75*pas;
                        xh = x_prec + pas * d;
                        [f_h,g_h,~, J_h] = rosenbrock(xh,params);
                    end
                end
                %sauvegarde des résultats intermédiares
                result.crit = [result.crit, f_h]; %critère par itération 
                result.grad = [result.grad, norm(g_h)]; %gradient par itération
                xval = [xval, xh]; %sauvegarde des itérées                

        end
        k=k+1;
    end
    %mise à jour de la variable result.stop
    if k == options.maxiter
        result.stop = 'Maxiter';
    elseif norm(g_h) < options.tolG
        result.stop = 'TolG';
    elseif norm((xh-x_prec)/xh) > options.tolX
        result.stop = 'TolX';
    elseif norm((f_h-f_prec)/f_h) > options.tolF
        result.stop = 'TolF';
    end
    result.iter = k;
    result.temp = toc;

end


