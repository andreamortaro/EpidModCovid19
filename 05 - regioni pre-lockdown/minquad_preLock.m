function L = minquad_preLock(K)     % funzionale da minimizzare

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

global  x0 tm ym Nass t_0 t_u tspan pnt

SI = @(t,x) [-K(1)*x(1)*x(2); K(1)*x(1)*x(2) - K(2)*x(2)]; % riscrivo sys con K

% Jac = @(t,x) [-K(1)*x(2), -K(1)*x(2);
%                K(1)*x(2),  K(1)*x(1) -  K(2)];
% options.Jacobian = Jac;
% [t, xm]  = eulerorosenbrock(SI,tspan,x0,options);

% difficile da controllare i tempi della soluzione
options.InitialStep = 1/pnt;    % così t(pnt*tm+1) deve coincidere con tm
[t, xm]  = rk4(SI,[t_0 t_u],x0,options);

if tm' ~= t(pnt*(tm-t_0)+1)     % devono coincidere
    warning("stai minimizzando nei punti sbagliati")
end

xm(:,3) = ones(length(t),1) - xm(:,1) - xm(:,2);        % ricavo R per post-processing
xm = Nass.*xm;                                          % ri-normalizzo da percentuale a Nass

% Calcolo numericamente la funzione dei minimi quadrati L
psi = 0.5;
phi = 1-psi;
n = 2;
deltat = t_u - t_0;
tt = pnt*tm+1;          % tempi soluzione (calcolata con diverso nstep da tm, piu punti)

L = 0;
for j = 2:length(tm)
    L = L+(psi*(ym(j,1)-xm(tt(j),2)).^n + phi*(ym(j,2)-xm(tt(j),3)).^n);  % misura minimi quadrati
end

L = deltat*L;

end
