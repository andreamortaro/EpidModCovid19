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

global  tstart tfinal tm ym Nass pnt

I = @(t) ym(1,1)*K(1)*exp(K(2)*t)./(K(1)+ym(1,1)*(exp(K(2)*t)-1));

% Calcolo numericamente la funzione dei minimi quadrati L
%phi = 0.985;
%psi = 1-phi;
%n = 2;
%deltat = tfinal-tstart;
%tt = pnt*(tm-tstart)+1;          % tempi soluzione (calcolata con diverso nstep da tm, piu punti)

L = 0;
for jj = 1:length(tm)
%     c = 100;
%     t0 = 9;    % tempo inizio
%     pot = 1./(1+exp(-c*(tm(jj)-t0)));
%     A = 100;
%     pot = A*pot;
    %pot = 1;
%     L = L+(phi*pot*(ym(jj,1)-xm(tt(jj),2)).^n + ...
%            psi*pot*(ym(jj,2)-xm(tt(jj),3)).^n);  % misura minimi quadrati
    L = L + (I(jj)-ym(jj,1)).^2;
end

mu = 0.001;
L = L/Nass + mu*(K(1)^0.5);


end
