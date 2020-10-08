function L = minquad_kdiscreti(K)     % funzionale da minimizzare

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

global  x0 tm ym Nass tl tr tspan pnt beta gamma mu

SEI = @(t,x) [-(beta - x(1)*x(3)/K)*x(1)*x(3);...
               (beta - x(1)*x(3)/K)*x(1)*x(3) - mu*x(2);...
               mu*x(2) - gamma*x(3)];
          
Jac = @(t,x) [ -beta*x(3) + 2*x(1)*(x(3)^2)/K, 0, -beta*x(1) + 2*(x(1)^2)*x(3)/K;
                beta*x(3) - 2*x(1)*(x(3)^2)/K, -mu,  beta*x(1) - 2*(x(1)^2)*x(3)/K;
                0,mu,-gamma];
options.Jacobian = Jac;

% la sol dipende da K
[t, xm]  = eulerorosenbrock(SEI,tspan,x0,options);
%options = odeset('Jacobian',Jac);
%[t, xm]  = ode23s(SEI,tspan,x0,options);


if tm' ~= t(pnt*(tm-tl)+1)
    warning("stai minimizzando nei punti sbagliati")
end

xm(:,4) = ones(length(t),1) - xm(:,1) - xm(:,2) - xm(:,3);	% ricavo R per post-processing
xm = Nass.*xm;                                      % ri-normalizzo da percentuale a Nass

% Calcolo numericamente la funzione dei minimi quadrati L
phi = 0.6;
psi = 0.25;
nu = 1-phi-psi;
n = 2;
tt = pnt*(tm-tl)+1;          % tempi soluzione (calcolata con diverso nstep da tm, piu punti)

L = 0;
for jj = 2:length(tm)
    L = L+(phi*(ym(jj,1)+ym(jj,2)-xm(tt(jj),2)-xm(tt(jj),3)).^n + ...
            psi*(ym(jj,1)-xm(tt(jj),2)).^n + ...
            nu*(ym(jj,3)-xm(tt(jj),4)).^n);  % misura minimi quadrati
end

%L = L/Nass;

end
