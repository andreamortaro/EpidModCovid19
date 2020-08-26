function [Kval,Lval] = plot_minquad_kdiscreti(kspan,pnt)     % funzionale da minimizzare

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

global t_u t_c Nass beta gamma Ibar Rbar

tm = t_u:1:t_c;
ym = [Ibar(tm+1),Rbar(tm+1)];

I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;          % dato iniziale in percentuale    

nstep = pnt*(t_c-t_u)+1;              % scritto cosi' per sapere dove minimizzare
tspan = linspace(t_u,t_c,nstep);      % intervallo dove minimizzare minquad

Kval = zeros(length(kspan),1);
Lval = zeros(length(kspan),1);

it = 0;
for K = kspan

    it=it+1;
    
    % simulo il modello con tale K
    SI = @(t,x) [-(beta - x(1)*x(2)/K)*x(1)*x(2);
                  (beta - x(1)*x(2)/K)*x(1)*x(2) - gamma*x(2)];
          
    Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K, -beta*x(1) + 2*(x(1)^2)*x(2)/K;
                    beta*x(2) - 2*x(1)*(x(2)^2)/K,  beta*x(1) - 2*(x(1)^2)*x(2)/K - gamma];
    options.Jacobian = Jac;
    
    [t, xm]  = eulerorosenbrock(SI,tspan,x0,options);

    if tm' ~= t(pnt*(tm-t_u)+1)
        warning("stai minimizzando nei punti sbagliati")
    end

    xm(:,3) = ones(length(t),1) - xm(:,1) - xm(:,2);	% ricavo R per post-processing
    xm = Nass.*xm;                                      % ri-normalizzo da percentuale a Nass

    
    % Calcolo numericamente la funzione dei minimi quadrati L
    psi = 0.5;
    phi = 1-psi;
    n = 2;
    deltat = t_c-t_u;
    tt = pnt*(tm-t_u)+1;          % tempi soluzione (calcolata con diverso nstep da tm, piu punti)

    L = 0;
    for j = 2:length(tm)
        L = L+(psi*(ym(j,1)-xm(tt(j),2)).^n + ...
               phi*(ym(j,2)-xm(tt(j),3)).^n);  % misura minimi quadrati
    end
    L = deltat*L;

    % salvo i dati
    Kval(it) = K;
    Lval(it) = L;
end


end
