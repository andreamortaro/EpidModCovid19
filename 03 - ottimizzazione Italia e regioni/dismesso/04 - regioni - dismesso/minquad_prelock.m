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

SI = @(t,x) [-K(1)*x(1)*x(2); K(1)*x(1)*x(2) - K(2)*x(2)]; % riscrivo sys con K

Jac = @(t,x) [-K(1)*x(2), -K(1)*x(2);
               K(1)*x(2),  K(1)*x(1) -  K(2)];
options.Jacobian = Jac;

nstep = pnt*(tfinal-tstart)+1;
tspan = linspace(tstart,tfinal,nstep);
[t, xm]  = eulerorosenbrock(SI,tspan,x0,options);

if tm' ~= t(pnt*(tm-tstart)+1)     % devono coincidere
    warning("stai minimizzando nei punti sbagliati")
end

xm(:,3) = ones(length(t),1) - xm(:,1) - xm(:,2);        % ricavo R per post-processing
xm = Nass.*xm;                                          % ri-normalizzo da percentuale a Nass

% Calcolo numericamente la funzione dei minimi quadrati L
phi = 0.5;
psi = 1-phi;
n = 2;
deltat = tfinal-tstart;
tt = pnt*(tm-tstart)+1;          % tempi soluzione (calcolata con diverso nstep da tm, piu punti)

L = 0;
for jj = 2:length(tm)
    L = L+(phi*(ym(jj,1)-xm(tt(jj),2)).^n + ...
           psi*(ym(jj,2)-xm(tt(jj),3)).^n);  % misura minimi quadrati
end

L = deltat*L;

end
