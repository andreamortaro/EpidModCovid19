function [days, K_disc]=stima_kdiscreti(data,kspan,window,K0_disc,pnt)

%
%   [days, K_disc] = stima_kdiscreti(data,kspan,window,K0_disc,options)
%
%   Per ogni punto di kspan integro su una finestra, dove i parametri sono
%   descritti nella struttura window, e ricavo il valore del parametro di controllo k.
%   Trovo un k per ogni punto (=giorno) di kspan.
%
%   INPUTS:
%   data        : struttura contenente i dati utili
%   kspan       : Intervallo temporale considerato.
%   window      : struttura che contiene i seguenti campi:
%                 - window.kl: estremo sinistro finestra integrazione rispetto t_i
%                 - window.kr: estremo destro finestra integrazione rispetto t_i
%                 - window.k0: guess iniziale
%   K0_disc     : guess iniziale
%
%   OUTPUTS:
%   days        : il giorno days(i) corrisponde al K_disc(i) parametro trovato
%   K_disc      : parametri discreti trovati per ogni punto di kpsan
%
%

global t_u t_c tl tr tm tspan ym Nass Ibar Rbar Ebar beta gamma mu x0 pnt  %#ok<REDEFGI>

% recupero i valori che servono
[Nass,Ibar,Rbar,~,Ebar] = data.value;
[~,~,t_u,t_c,~] = data.time;
[beta,gamma,mu] = data.parameters;

h = window.h;
kl = window.kl;
kr = window.kr;

% aggiungo nodi per migliore risoluzione minquad
if (nargin == 3)
  pnt = 1;
end

days    = zeros(t_c-kr*h-t_u,1);
K_disc  = zeros(t_c-kr*h-t_u,1);

it = 1;

if iscolumn(kspan)
    kspan = kspan';
end

for t_i=kspan
    
    tl = t_i-kl*h;                  % tempo iniziale
    tr = t_i+kr*h;                  % tempo finale
    nstep = pnt*(tr-tl)+1;          % scritto cosi' per sapere dove minimizzare
    tspan = linspace(tl,tr,nstep);  % intervallo dove minimizzare minquad
    
    tm = tl:1:tr;
    ym = [Ebar(tm+1),Ibar(tm+1),Rbar(tm+1)];
    
    E0 = Ebar(tl+1); I0 = Ibar(tl+1); R0 = Rbar(tl+1); S0 = Nass-E0-I0-R0;
    x0 = [S0;E0;I0]/Nass;          % dato iniziale in percentuale

    problem.options = optimoptions('fmincon','Display','iter');
    problem.solver = 'fmincon';
    problem.objective = @minquad_kdiscreti;     % funzionale obiettivo minimizzare
    problem.x0 = K0_disc;                       % guess iniziale
    problem.lb = 0;                             % lower bound
    problem.nonlcon = @(K)mycon(K);             % vincolo non lineare su k (=beta>0)

    % salvo giorno e parametro ottenuto
    days(it) = t_i;
    K_disc(it) = fmincon(problem);     % salvo il parametro ottenuto
    
    K0_disc = K_disc(it);                 % update guess iniziale
    it = it+1;                      % update iterazioni
end

end


% VINCOLO NON LINEARE SU K

function [c,ceq] = mycon(K)

global x0 beta gamma mu tl tr

    SEI = @(t,x) [-(beta - x(1)*x(3)/K)*x(1)*x(3);...
                   (beta - x(1)*x(3)/K)*x(1)*x(3) - mu*x(2);...
                    mu*x(2)- gamma*x(3)];
          
    Jac = @(t,x) [ -beta*x(3) + 2*x(1)*(x(3)^2)/K, 0, -beta*x(1) + 2*(x(1)^2)*x(3)/K;...
                    beta*x(3) - 2*x(1)*(x(3)^2)/K, -mu,  beta*x(1) - 2*(x(1)^2)*x(3)/K;...
                    0, mu,-gamma];
    options.Jacobian = Jac;    
    
    nstep = (tr-tl)+1;
    [~, xm]  = eulerorosenbrock(SEI,linspace(tl,tr,nstep),x0,options);
%    [~, xm]  = ode15s(SEI,linspace(tl,tr,nstep),x0);


    % controllo la condizione con i valori in percentuale, altrimenti e'
    % impossibile verificare la condizione
    
	% Nonlinear inequality constraints (c(K)<=0)
    c = xm(:,1).*xm(:,3) - beta*K;   % deve essere <=0 (ATTENZIONE all'=)
    
    ceq = [];
end
