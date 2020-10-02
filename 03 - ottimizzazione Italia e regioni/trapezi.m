function [tout, yout] = trapezi(odefun, tspan, y0, options)

%
%   [tout, yout] = trapezi (odefun, tspan, y0, options)
%
%   Metodo dei trapezi
%
%   Risolve l'ODE y'(t) = odefun(t,y(t))
%
%   INPUTS:
%   odefun  : campo vettoriale associato
%   tspan   : vettore di tempi in cui vogliamo calcolare la soluzion
%   y0      : dato iniziale
%   options : variabile di tipo struttura:  - options.Jacobian = d(odefun)/dy
%
%   OUTPUTS:
%   tout    : vettore di tempi nei quali abbiamo calcolato la soluzione
%   yout    : matrice le cui colonne sono le soluzioni ad ogni tempo
%

m = length(tspan)-1;                % numero di passi
yout = zeros(length(y0), m+1);
yout(:,1) = y0(:);                  % y0 in prima colonna yout
k = tspan(2)-tspan(1);              % calcolo il passo (costante)
dodefundy = options.Jacobian;

tol = k^2/100; % tolleranza proporzionale all'errore che pensiamo di commettere

for n = 1:m
    
    F = @(x) x-yout(:,n)-k/2*odefun(tspan(n),yout(:,n))-k/2*odefun(tspan(n+1),x);
    J = @(x) eye(length(x))-k/2*dodefundy(tspan(n+1),x);
    yout(:,n+1)=yout(:,n);                        % guess iniziale per il metodo di Newton
    delta = -J(yout(:,n+1)) \ F(yout(:,n+1)) ;    % incognita del metodo di Newton
    
    % Metodo di Newton
    while (norm(delta) > tol)
        yout(:,n+1) = yout(:,n+1)+delta;
        delta = -J(yout(:,n+1)) \ F(yout(:,n+1)) ;
    end
    
    yout(:,n+1) = yout(:,n+1)+delta;     % l'ultima correzione che rende la soluzione più precisa

end

tout = tspan(:);
yout = yout.';

end
