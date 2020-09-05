function [days, k_c]=stima_kdiscreti(kspan,window,k0_c,options)

%
%   [iter, days, k_c] = stima_kdiscreti(kspan,window,k0_c,options)
%
%   Per ogni punto di kspan integro su una finestra, dove i parametri sono
%   descritti nella struttura window, e ricavo il valore del parametro di controllo k.
%   Trovo un k per ogni punto (=giorno) di kspan.
%
%   INPUTS:
%   kspan       : Intervallo temporale considerato.
%   window      : struttura che contiene i seguenti campi:
%                 - window.kl: estremo sinistro finestra integrazione rispetto t_i
%                 - window.kr: estremo destro finestra integrazione rispetto t_i
%                 - window.k0: guess iniziale
%
%   OUTPUTS:
%   iter       : numero di iterazioni fatte nel ciclo for
%                (corrisponde con length(kspan) - kr*h)
%   days       : il giorno days(i) corrisponde al k_c(i) parametro trovato
%   k_c        : parametri discreti trovati per ogni punto di kpsan
%
%


global x0 tm ym Nass tl tr tspan pnt t_c t_u Ibar Rbar

h = window.h;
kl = window.kl;
kr = window.kr;

pnt = 1;        % default
if (nargin == 4) 
  if (isfield(options,'pnt'))   % aggiungo nodi per migliore risoluzione minquad
    pnt = options.pnt;
  end
end

days = zeros(t_c-kr*h-t_u,1);
k_c  = zeros(t_c-kr*h-t_u,1);

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
    ym = [Ibar(tm+1),Rbar(tm+1)];
    
    I0 = Ibar(tl+1); R0 = Rbar(tl+1); S0 = Nass-I0-R0;
    x0 = [S0;I0]/Nass;          % dato iniziale in percentuale
                                % x0 serve per fmincon(problem), dentro risolvo una ode
    
    % Minimizzazione: nella finestrella di 7 giorni centrata in t_i cerchiamo
    % il parametro k_c di controllo

    %problem.options = optimoptions('fmincon','Display','iter',...
    %                   'PlotFcn',{@optimplotfval,@optimplotfirstorderopt});
    problem.options = optimoptions('fmincon','Display','iter');
    problem.solver = 'fmincon';
    problem.objective = @minquad_kdiscreti;     % funzionale obiettivo minimizzare
    problem.x0 = k0_c;                          % guess iniziale min
    problem.lb = 0;                             % lower bound
    problem.nonlcon = @(K)mycon(K);             % vincolo non lineare su k (=beta>0)
    
    %prove
    problem.OptimalityTolerance = 1e-12;
    problem.StepTolerance = 1e-14;
    problem.FunctionTolerance = 1e-18;
    problem.ConstraintTolerance = 1e-12;
    %problem.Algorithm = 'active-set';
    problem.MaxFunctionEvaluations = 5000; %max per fmincon con punto interno 3000
    problem.MaxIterations = 1500;
    
    
    % salvo giorno e parametro ottenuto
    days(it) = t_i;
    k_c(it) = fmincon(problem);     % salvo il parametro ottenuto
    
    k0_c = k_c(it);                 % update guess iniziale
    it = it+1;                      % update iterazioni
end

end


% VINCOLO NON LINEARE SU K

function [c,ceq] = mycon(K)

global x0 beta gamma tspan Nass tl tr
    
    SI = @(t,x) [-(beta - x(1)*x(2)/K)*x(1)*x(2);
                  (beta - x(1)*x(2)/K)*x(1)*x(2) - gamma*x(2)];
          
    Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K, -beta*x(1) + 2*(x(1)^2)*x(2)/K;
                    beta*x(2) - 2*x(1)*(x(2)^2)/K,  beta*x(1) - 2*(x(1)^2)*x(2)/K - gamma];
    options.Jacobian = Jac;
    
    tt = linspace(tl,tr,11);
    [t, xm]  = eulerorosenbrock(SI,tt,x0,options);
    
    %[t, xm]  = eulerorosenbrock(SI,tspan,x0,options);
    
    %xm(:,3) = ones(length(t),1) - xm(:,1) - xm(:,2);      % ricavo R per post-processing
    %xm = Nass.*xm;                                       % ???
    
    % controllo la condizione con i valori in percentuale, altrimenti e'
    % impossibile verificare la condizione
    
	% Nonlinear inequality constraints (c(K)<=0)
    c = xm(:,1).*xm(:,2) - beta*K;   % deve essere <=0 (ATTENZIONE all'=)
    
    ceq = [];
end
