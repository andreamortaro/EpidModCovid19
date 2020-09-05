function [tout,yout] = rk4(odefun,tspan,y0,options)

%
%   [tout,yout] = rk4(odefun,tspan,y0,options)
%
%   Metodo di Eulero modificato (esplicito) (Runge-Kutta a 2 stadi)
%
%   Risolve l'ODE y'(t) = odefun(t,y(t))
%
%   INPUTS:
%   odefun  : campo vettoriale associato
%   tspan   : vettore di tempi [t0,t0+tstar]
%   y0      : dato iniziale
%   options : variabile di tipo struttura:  - options.InitialStep   = passo
%
%   OUTPUTS:
%   tout    : vettore di tempi nei quali abbiamo calcolato la soluzione
%   yout    : matrice le cui colonne sono le soluzioni ad ogni tempo
%
%   Il passo iniziale InitialStep viene posto di default a 0.5 se non è
%   presente nel campo options
%

InitialStep = 0.05;
if (nargin == 4) 
  if (isfield(options,'InitialStep'))
    InitialStep = options.InitialStep;
  end
end
k = InitialStep;

% Tableau
A(2,1) = 1/2;
A(3,2) = 1/2;
A(4,3) = 1;
c(2) = 1/2;
c(3) = 1/2;
c(4) = 1;
b(1) = 1/6;
b(2) = 1/3;
b(3) = 1/3;
b(4) = 1/6;

v = 4;  % numero di stadi

n = 1;
tout(n) = tspan(1);     % tempo iniziale
yout(:,n) = y0(:);      % stima iniziale

xi = zeros(length(y0),v);
f = zeros(length(y0),v);

while (tspan(2)-tout(n) > eps)      % finche non raggiungo t finale a meno di un errore macchina eps
    k = min(k,tspan(2)-tout(n));    % controllo la distanza dalla fine e diminuisco k se necessario

    % calcolo xi e f
    for i=1:v
        xi(:,i) = yout(:,n);
        for j = 1:i-1
            xi(:,i) = xi(:,i)+k*A(i,j)*f(:,j);
        end
        f(:,i)=odefun(tout(n)+c(i)*k, xi(:,i));
    end

    % calcolo la soluzione
    yout(:,n+1) = yout(:,n);
    for j=1:v
        yout(:,n+1) = yout(:,n+1)+k*b(j)*f(:,j);
    end

    tout(n+1) = tout(n)+k;
    n = n+1;
end

% stranezza:

tout = tout.';
yout = yout.';

% se l'ultimo punto è troppo ravvicinato al precedente lo taglio
if(tout(end)-tout(end-1)<InitialStep/100)
    tout = tout(1:end-1);
    yout = yout(1:end-1,:);
end

