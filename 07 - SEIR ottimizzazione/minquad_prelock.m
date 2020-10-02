function L = minquad_prelock(K)     % funzionale da minimizzare

%
%   L = minquad(K)
%
%   Definisco il funzionale L da minimizzare
%
%   INPUTS:
%   K       : incognita
%
%   OUTPUTS:
%   L       : funzione minimi quadrati
%

global  tstart tfinal tm ym x0 Nass pnt

SEI = @(t,x) [-K(1)*x(1)*x(3);...
              K(1)*x(1)*x(3) - K(3)*x(2)
              K(3)*x(2)- K(2)*x(3)]; % riscrivo sys con K

Jac = @(t,x) [-K(1)*x(3),0, -K(1)*x(3);...
               K(1)*x(3), -K(3), K(1)*x(1);...
               0,K(3),-K(2)];
options.Jacobian = Jac;

nstep = pnt*(tfinal-tstart)+1;
tspan = linspace(tstart,tfinal,nstep);
[t, xm]  = eulerorosenbrock(SEI,tspan,x0,options);

if tm' ~= t(pnt*(tm-tstart)+1)     % devono coincidere
    warning("stai minimizzando nei punti sbagliati")
end

xm(:,4) = ones(length(t),1) - xm(:,1) - xm(:,2) - xm(:,3);        % ricavo R per post-processing
xm = Nass.*xm;                                          % ri-normalizzo da percentuale a Nass

% Calcolo numericamente la funzione dei minimi quadrati L
phi = 0.6;
psi = 0.3;
nu = 1-phi-psi;
n = 2;
deltat = tfinal-tstart;
tt = pnt*(tm-tstart)+1;          % tempi soluzione (calcolata con diverso nstep da tm, piu punti)

L = 0;
for jj = 2:length(tm)
    c = 100;
    t0 = 9;    % tempo inizio
    pot = 1-1./(1+exp(c*(tm(jj)-t0)));
    A = 5;
    pot = A*pot;
    %pot = 1;
    L = L+(phi*pot*(ym(jj,1)-xm(tt(jj),2)).^n + ...
           psi*pot*(ym(jj,2)-xm(tt(jj),3)).^n + ...
            nu*pot*(ym(jj,3)-xm(tt(jj),4)).^n);  % misura minimi quadrati
end

mu = 1;
L = L/Nass + mu*(1*K(1)^2 + 0*K(2)^2);

%L = L/Nass + mu*(K(1)^0.5);


end
