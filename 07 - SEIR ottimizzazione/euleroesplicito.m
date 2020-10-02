function [tout, yout] = euleroesplicito (odefun, tspan, y0)

%
%   [tout, yout] = euleroesplicito (odefun, tspan, y0)
%
%   Risolve l'ODE y'(t) = odefun(t,y(t))
%
%   INPUTS:
%   odefun  : campo vettoriale associato
%   tspan   : vettore di tempi in cui vogliamo calcolare la soluzione
%   y0      : dato iniziale
%
%   OUTPUTS:
%   tout    : vettore di tempi nei quali abbiamo calcolato la soluzione
%   yout    : matrice le cui colonne sono le soluzioni ad ogni tempo
%

m = length(tspan) - 1;
yout = zeros(length(y0),m+1);      % inizializzo perche' costoso allocare vettori ogni volta
yout(:,1) = y0(:);                 % y0 nella prima colonna di yout
k = tspan(2) - tspan(1);           % calcolo il time-step k

% se tspan non e' equispaziato posso cambiare k ad ogni iterazione

for n = 1:m     % parto dalla seconda colonna perche' ho y0 nella prima
	yout(:, n + 1) = yout(:,n) + k * odefun(tspan(n), yout(:,n));
end

tout = tspan(:);    % tout colonna
yout = yout.';      % trasposizione reale

end
