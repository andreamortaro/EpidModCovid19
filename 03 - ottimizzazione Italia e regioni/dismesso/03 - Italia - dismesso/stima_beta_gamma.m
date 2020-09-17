function [t,x,beta,gamma] = stima_beta_gamma(data,K0,pnt)

%
%   [t, x, beta, gamma] = prelock(data, K0, options)
%
%   Stima di beta e gamma nell'intervallo [data(1).prelock,data(2).prelock]
%   e simulazione del modello
%
%   INPUTS:
%   data        : struttura che contiene i dati utili come Ibar e Rbar
%   K0          : guess ottimizzazione per trovare beta-gamma.
%
%   OUTPUTS:
%   t           : tempi soluzione del modello simulato
%   x           : soluzione modello simulato
%   beta        : parametro del SIR ottenuto dall'ottimizzazione
%   gamma       : parametro del SIR ottenuto dall'ottimizzazione
%

global tstart tfinal date Ibar Rbar Nass x0 tm ym tspan pnt  %#ok<REDEFGI>

% recupero i valori che servono

[tstart,tfinal] = data.prelock;

[Nass,Ibar,Rbar] = data.value;
date = data(5).time;

tm  = tstart:1:tfinal;
ym = [Ibar(tm+1),Rbar(tm+1)];

nstep = pnt*(tfinal-tstart)+1;
tspan = linspace(tstart,tfinal,nstep);

I0 = Ibar(tstart+1); R0 = Rbar(tstart+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;                           % dato iniziale in percentuale

% Minimizzazione: cerco beta e gamma

problem.options     = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
problem.solver      = 'fmincon';
problem.objective   = @minquad_prelock;
problem.x0 = K0;
problem.lb = [0,0];                 % impongo positivi i parametri cercati
problem.OptimalityTolerance = 1e-12;
problem.StepTolerance = 1e-12;
problem.FunctionTolerance = 1e-12;

K = fmincon(problem)

beta = K(1); gamma = K(2);  %  parametri stimati

%% update sys, nuovi parametri stimati
SI = @(t,x) [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma*x(2)];

Jac = @(t,x) [-beta*x(2), -beta*x(2);
               beta*x(2), beta*x(1) - gamma];
opt.Jacobian = Jac;
opt.InitialStep = 0.001;

[t, x]  = rk4(SI,[tstart,tfinal],x0,opt);             % simulazione modello

x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R per post-processing
x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass

end

