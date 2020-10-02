function xsol = odesol(fun,t,x)

%
%   [xsol = odesol(fun,t,x)
%
%   Risolve l'ODE y'(t) = fun(t)
%
%   INPUTS:
%   fun         : campo vettoriale associato
%   t           : vettore di tempi in cui vogliamo calcolare la soluzione
%   x           : dato iniziale
%
%   OUTPUTS:
%   xsol        : matrice le cui colonne sono le soluzioni ad ogni tempo
%
%   NOTA: per analisi di sensitivita'
%

M = length(t);
xsol(1,:) = x;  % salvo il guess

%Jn = options.Jacobian;

for m=1:M-1
    h = t(m+1)-t(m);            % step temporale
    %x1 = x + h*phi1m(h*Jn(m,x))*fun(m, x);          % eulero-rosenbrock
    x1 = x + h*fun(m,x);       % stima successiva
    x  = x1;
    xsol(m+1,:) = x;            % inserisco la soluzione nei risultati
end

end
