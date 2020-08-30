function A = fitting_k(A0,days,k_c)

%
%   L = minquad_kcontinuo(A)
%
%   poi per t~days dovra valere K(t) = k_c(days)
%
%   INPUTS:
%   A0      : guess iniziale per i parametri per il fitting
%   days    : giorni dei parametri da fittare
%   k_c     : parametri da fittare
%
%   OUTPUTS:
%   A       : incognita della forma A = [A(1),A(2),A(3)]
%

a = A0(1); b = A0(2); c = A0(3);

problem2.options    = optimoptions('fmincon','Display','iter');
problem2.solver     = 'fmincon';
problem2.objective  = @minquad_kcontinuo;           % funzionale obiettivo minimizzare
problem2.x0         = [a,b,c];                      % guess iniziale
%problem2.nonlcon = @(A)mycon_discreto(A);                  % vincolo non lineare su k (=beta>0)

% tecnicamente dovrei imporre il vinc non lin, tuttavia
%l'ho imposto quando ho cercato i k discreti ???

A = fmincon(problem2);

end

% VINCOLO NON LINEARE SU K Quando voglio fittarei k discreti con una
% function k(t) continua

function [c,ceq] = mycon_discreto(A)

global x0 beta gamma t_u t_c Ibar Rbar Nass

    K = @(t) A(1)*exp(A(2)*t).*(1-exp(A(3)*t));
 
    SI = @(t,x) [-(beta - x(1)*x(2)/K(t))*x(1)*x(2);
                  (beta - x(1)*x(2)/K(t))*x(1)*x(2) - gamma*x(2)];
          
    Jac = @(t,x) [ -beta*x(2) + 2*x(1)*(x(2)^2)/K(t), -beta*x(1) + 2*(x(1)^2)*x(2)/K(t);
                    beta*x(2) - 2*x(1)*(x(2)^2)/K(t),  beta*x(1) - 2*(x(1)^2)*x(2)/K(t) - gamma];
    options.Jacobian = Jac;
    
    I0 = Ibar(t_u+1); R0 = Rbar(t_u+1); S0 = Nass-I0-R0;
    x0 = [S0;I0]/Nass;          % dato iniziale in percentuale
    
    tspan = linspace(t_u,t_c,74);            % cosi non ho problemi con matrici singolari
    [t, xm]  = eulerorosenbrock(SI,tspan,x0,options);

    %xm = Nass.*xm;                                        % ???
    
    % controllo con valori in percentuale, altrimenti impossibile verificare la condizione
    
	% Nonlinear inequality constraints (c(K)<=0)
    
    if t ~= tspan'
        warning('attenzione ai vincoli')
    end
    
    c = xm(:,1).*xm(:,2) - 2*beta*K(t);   % deve essere <=0 (ATTENZIONE all'=)
    ceq = [];
end