close all
clear all

% DEFINITION DES PARAMETRES DU PROBLEME
params.fonction='rosenbrock';


% DEFINITION DE L ALGORITHME D OPTIMISATION 
%options.method='gradient';
% options.pas='variable';
options.tolG= 1e-8;
options.maxiter = 1e3;
options.tolX= 1e-8;
options.tolF = 1e-8;


x0 = [-4;10];
% --------------------------affichage de g----------------------
options.method='gradient';
options.pas='variable';
[xh_g,result_g,~] = optimdescent(params.fonction,2,options,x0);
figure(1)
plot(result_g.grad,'r-o')
xlabel('itération')
ylabel('norme du gradient de f')

%affichage du gradient conjuge
options.pas='variable';
options.method='gradient conjuge';
[xh_gc,result_gc,~] = optimdescent(params.fonction,2,options,x0); 
figure(1)
hold on 
plot(result_gc.grad,'b-x')

%affichage du Newton pas unitaire
options.pas='fixe';
options.method='Newton';
[xh_Nf,result_Nf,~] = optimdescent(params.fonction,2,options,x0);
hold on 
plot(result_Nf.grad,'g-s')

%affichage du Newton pas variable
options.pas='variable';
options.method='Newton';
[xh_Nv,result_Nv,~] = optimdescent(params.fonction,2,options,x0);
hold on 
plot(result_Nv.grad,'c-v')

%affichage du Quasi-Newton
options.pas='variable';
options.method='Quasi-Newton';
[xh_qn,result_qn,~] = optimdescent(params.fonction,2,options,x0);
hold on 
plot(result_qn.grad,'k-')
legend('gradient', 'gradient conjuge', 'Newton pas unitaire', 'Newton pas variable', 'Quasi-Newton')

% -----------------affichage f-----------------------------
options.method='gradient';
options.pas='variable';
[xh_g,result_g,~] = optimdescent(params.fonction,2,options,x0);
figure()
plot(result_g.crit,'r-o')
xlabel('itération')
ylabel('critère de f')

%affichage du gradient conjuge
options.pas='variable';
options.method='gradient conjuge';
[xh_gc,result_gc,~] = optimdescent(params.fonction,2,options,x0);
hold on 
plot(result_gc.crit,'b-x')

%affichage du Newton pas unitaire
options.pas='fixe';
options.method='Newton';
[xh_Nf,result_Nf,~] = optimdescent(params.fonction,2,options,x0);
hold on 
plot(result_Nf.crit,'g-s')

%affichage du Newton pas variable
options.pas='variable';
options.method='Newton';
[xh_Nv,result_Nv,~] = optimdescent(params.fonction,2,options,x0);
hold on 
plot(result_Nv.crit,'c-v')

%affichage du Quasi-Newton
options.pas='variable';
options.method='Quasi-Newton';
[xh_qn,result_qn,~] = optimdescent(params.fonction,2,options,x0);
hold on 
plot(result_qn.crit,'k-')
legend('gradient', 'gradient conjuge', 'Newton pas unitaire', 'Newton pas variable', 'Quasi-Newton')

% 
% 
% 
% 
% 
% 
% 
% % AFFICHAGE DES RESULTATS
% x1 = -4:0.05:4;
% x2 = -5:0.05:20;
% %initialisation du vecteur f_0 pour le tracé
% f_0 = zeros(length(x2),length(x1));
% %calcul de la fonction objectif en chaque point de l'intervalle x1 et x2
% for i=1:length(x1)
%     for j=1:length(x2)
%         x = [x1(i) ; x2(j)];
%         [f, ~, ~] = rosenbrock(x,2);
%         f_0(j,i) = f;
%     end
% end
% %grillage 2D de points sur lesquels la fonction f_0 sera évaluée
% [X1, X2] = meshgrid(x1, x2);
% %tracé des lignes de niveaux de la fonction f_0 sur une autre figure 
% figure();
% contour(X1,X2,f_0,100)
% 
% xlabel('x1')
% ylabel('x2')
% title('lignes de niveaux de la fonction de Rosenbrock pour n=2 et b=2')
% colorbar
% 
% %affichage de la trajectoire des itérées
% hold on
% plot(xval(1,:),xval(2,:),'k-x')
% 
% 
% %affichage de la norme du gradient à chaque itération
% % 
% % figure()
% % plot(result.grad,'r-o')
% % xlabel('itération')
% % ylabel('norme du gradient de f')
% % figure()
% % plot(result.crit,'r-o')








% A COMPLETER
