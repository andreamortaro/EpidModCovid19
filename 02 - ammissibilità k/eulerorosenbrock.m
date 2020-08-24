function [tout, yout] = eulerorosenbrock (odefun, tspan, y0, options)

%
%   [tout, yout] = eulerorosenbrock (odefun, tspan, y0, options)
%
%   Risolve l'ODE y'(t) = odefun(y(t)) (sistema autonomo)
%
%   INPUTS:
%   odefun      : campo vettoriale associato
%   tspan       : vettore di tempi in cui vogliamo calcolare la soluzione
%   y0          : dato iniziale
%   options     : variabile di tipo struttura:  - options.Jacobian = d(odefun)/dy
%
%   OUTPUTS:
%   tout    : vettore di tempi nei quali abbiamo calcolato la soluzione
%   yout    : matrice le cui colonne sono le soluzioni ad ogni tempo
%

m = length(tspan)-1;                    % numero di passi
yout = zeros(length(y0), m+1);
yout(:,1) = y0(:);                      % y0 in prima colonna yout
k = tspan(2)-tspan(1);                  % calcolo il passo (costante)

Jn = options.Jacobian;

for n = 1:m
    yout(:,n+1) = yout(:,n) + k*phi1m(k*Jn(tspan(n),yout(:,n)))*odefun(tspan(n), yout(:,n));
end

tout = tspan(:);
yout = yout.';

end
