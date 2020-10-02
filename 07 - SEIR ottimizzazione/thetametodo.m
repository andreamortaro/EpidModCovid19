function [tout, yout] = thetametodo (odefun, tspan, y0, options)

%
%   [tout, yout] = thetametodo (odefun, tspan, y0, options)
%
%   Risolve l'ODE y'(t) = odefun(t,y(t))
%
%   INPUTS:
%   odefun  : campo vettoriale associato
%   tspan   : vettore di tempi in cui vogliamo calcolare la soluzione
%   y0      : dato iniziale
%   options : variabile di tipo struttura: - options.Theta    = theta
%                                          - options.Jacobian = d(odefun)/dy
%   OUTPUTS:
%   tout    : vettore di tempi nei quali abbiamo calcolato la soluzione
%   yout    : matrice le cui colonne sono le soluzioni ad ogni tempo
%
%   Note:
%   theta = 0     --> Eulero esplicito
%   theta = 1/2   --> Trapezi
%   theta = 1     --> Eulero implicito
%

m = length(tspan)-1;                % numero di passi
yout = zeros(length(y0), m+1);
yout(:,1) = y0(:);                  % y0 nella prima colonna di yout
k = tspan(2) - tspan(1);              % calcolo il passo (costante)
theta = options.Theta;
dodefundy = options.Jacobian;

tol = k^2/100; % tolleranza proporzionale all'errore che pensiamo di commettere

for n = 1:m   
    F = @(x) x-yout(:,n)-k*(1-theta)*odefun(tspan(n),yout(:,n))- ...
        k*theta*odefun(tspan(n+1),x);
    J = @(x) eye(length(x))-k*theta*dodefundy(tspan(n+1),x);
    yout(:,n+1)=yout(:,n);                        % guess iniziale per il metodo di Newton
    delta = -J(yout(:,n+1)) \ F(yout(:,n+1)) ;    % incognita del metodo di Newton
    
    % Metodo di Newton
    while (norm(delta) > tol)
        yout(:,n+1) = yout(:,n+1)+delta;
        delta = -J(yout(:,n+1)) \ F(yout(:,n+1)) ;
    end
    
    yout(:,n+1) = yout(:,n+1)+delta;     % correzione che rende la soluzione più precisa

end

tout = tspan(:);
yout = yout.';

end

% se sbaglio di poco la derivata allora Newton converge ancora, ma e' piu' lento
% MEMENTO: per Newton il delta deve decadere con il quadrato ad ogni passo
