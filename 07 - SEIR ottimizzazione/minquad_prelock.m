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
%[t, xm]  = eulerorosenbrock(SEI,tspan,x0,options);
[t, xm]  = ode15s(SEI,tspan,x0);

if tm' ~= t(pnt*(tm-tstart)+1)     % devono coincidere
    warning("stai minimizzando nei punti sbagliati")
end

xm(:,4) = ones(length(t),1) - xm(:,1) - xm(:,2) - xm(:,3);        % ricavo R per post-processing
xm = Nass.*xm;                                          % ri-normalizzo da percentuale a Nass

% Calcolo numericamente la funzione dei minimi quadrati L
phi = 0.55;
psi = 0.35;
nu = 1-phi-psi;
n = 2;
deltat = tfinal-tstart;
tt = pnt*(tm-tstart)+1;          % tempi soluzione (calcolata con diverso nstep da tm, piu punti)

L = 0;
for jj = 2:length(tm)
    c = 100;
    t0 = 20;    % tempo inizio
    pot = 1./(1+exp(-c*(tm(jj)-t0)));
    A = 10;
    pot = A*pot;
    %pot = 1;
    L = L+(phi*pot*(ym(jj,1)+ym(jj,2)-xm(tt(jj),2)-xm(tt(jj),3)).^n + ...
            psi*pot*(ym(jj,1)-xm(tt(jj),2)).^n + ...
            nu*pot*(ym(jj,3)-xm(tt(jj),4)).^n);  % misura minimi quadrati
end

mmu = 1;
L = L/Nass + mmu*(1*K(1)^2 + 0*K(2)^2 + 1*K(3)^2);

%L = L/Nass + mmu*(K(1)^0.5);


end
