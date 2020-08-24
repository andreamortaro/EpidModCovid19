%{
Dinamica cambiata: il 09/03/20 parte il lockdown e nel sir entra in gioco
il fattore di controllo K.
Aggiustamento versione 2020-06-01: trovo un K per ogni giorno di lockdown e
ricavo un K(t) con t=giorno. Provo poi a fittare tale K(t) con una funzione
continua da inserire in un secondo momento nel modello.
%}

%% Identificazione parametri:

close all
clear all
clc

global x0 tm ym Nass tl tr tspan pnt beta gamma k_history

% Upload dati protezione civile
tmp = fullfile('..','..','..','0_dpc_data','dati-andamento-nazionale');
[status,result] = fileattrib(tmp);
path_folder = result.Name;              % percorso alla cartella
[date,Ibar,Rbar] = data_read_dpc(path_folder);
% Attenzione: poiche' gli indici in matrice partono da 1 e ho considerato
% l'istante iniziale pari a 0
% date(i), Ibar(i), Rbar(i) e' in corrispondenza con t_i-1

% DATI:
t_u = 14;                          % 2020-03-09: inizio lockdown
t_c = length(date)-1;              % decremento perche' parto da 0

% parametri trovati prima con main.m
beta    = 0.30883;
gamma = 0.04953;

Nass = 60317000;            % popolazione italiana istat 11.02.2020


h = 1; % daily time step
% finestre da 7 giorni
k_l = 3;
k_r = 4;
k0_c = 0.01;     % guess iniziale con cui far partire la prima minimizzazione

% variando guess:
% - in [0.05,1e5] ottengo gli stessi k_c
% - in [1e-8,0.01] per i primi 10-20 giorni ottengo dei k tendenti a 0
%
% mi sembrano piu convincenti i valori nel primo intervallo

% per salvare le prove con vari guess
txt = sprintf('guess=%g.txt',k0_c);
fileID = fopen(txt,'w');

fig = sprintf('fig_guess=%g.pdf',k0_c);

k_history = zeros(t_c-k_r*h-t_u,3);
it = 1;

for t_i=(t_u:1:t_c-k_r*h)
    
    it
    
    tl = t_i-k_l*h;     % tempo iniziale
    tr = t_i+k_r*h;     % tempo finale
    pnt = 20;            % aggiungo nodi per migliore risoluzione minquad
    nstep = pnt*(tr-tl)+1;      % scritto cosi' per sapere dove minimizzare
    tspan = linspace(tl,tr,nstep);
    
    tm = tl:1:tr;
    ym = [Ibar(tm+1),Rbar(tm+1)];
    
    I0 = Ibar(tl+1); R0 = Rbar(tl+1); S0 = Nass-I0-R0;
    x0 = [S0;I0]/Nass;          % dato iniziale in percentuale
    % x0 serve per fmincon(problem), dentro risolvo una ode
    
    % Minimizzazione: nella finestrella di 7 giorni centrata in t_i cerchiamo
    % il parametro k_c di controllo

    % problem.options = optimoptions('fmincon','Display','iter');
    problem.options = optimoptions('fmincon');
    problem.solver = 'fmincon';
    problem.objective = @minquad_window;       % funzionale obiettivo minimizzare
    problem.x0 = k0_c;                          % guess iniziale min
    problem.lb = 0;                             % lower bound
    problem.nonlcon = @(K)mycon(K);             % vincolo non lineare su k (=beta>0)
    
    %fmincon(problem)        % stampo il risultato dell'ottimizzazione
    
    % Update sistema e parametro strimato
    
    k_c = fmincon(problem);                 % update parametri stimati
    k_history(it,:) = [it, t_i, k_c];         % salvo il parametro ottenuto
    k0_c = k_history(it,3);                   % aggiorno il dato iniziale

    % sistema aggiornato
    SI = @(t,x) [-(beta - 0.5*x(1)*x(2)/k_c)*x(1)*x(2);
                  (beta - 0.5*x(1)*x(2)/k_c)*x(1)*x(2) - gamma*x(2)];

    Jac = @(t,x) [ -beta*x(2) + x(1)*(x(2)^2)/k_c, -beta*x(1) + (x(1)^2)*x(2)/k_c;
                    beta*x(2) - x(1)*(x(2)^2)/k_c,  beta*x(1) - (x(1)^2)*x(2)/k_c - gamma];
    options.Jacobian = Jac;

    [t, x]  = eulerorosenbrock(SI,tspan,x0,options);    % risolvo il sistema

    x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R post-processing
    x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass

    it = it+1;  % update iterazioni
end

T = table(k_history(:,1),k_history(:,2),k_history(:,3),'VariableNames',{'iter' 't_i' 'k_c(t_i)'}) %#ok<NOPTS>

M = [(1:1:length(k_history(:,1)));k_history(:,3)'];
fprintf(fileID,'%i      %i      %f\n', k_history');

figure(1)
plot(k_history(:,2),k_history(:,3),'*')
hold on

exportgraphics(figure(1),fig,'ContentType','vector',...
               'BackgroundColor','none')


% Minimizzazione: Fitto i k_c discreti ottenuti e ricavo k(t)

a = 0.1;
b = 0.1;
c = 0.1;

problem2.options = optimoptions('fmincon','Display','iter');
problem2.solver = 'fmincon';
problem2.objective = @minquad_k;             % funzionale obiettivo minimizzare
problem2.x0 = [a,b,c];                    % guess iniziale
%problem2.lb = 0;                             % lower bound
%problem2.nonlcon = @(K)mycon(K);             % vincolo non lineare su k (=beta>0)

A = fmincon(problem2)
a = A(1); b = A(2); c = A(3);

% update function
k_fun = @(t) a*exp(b*t).*(1-exp(c*t));

tt = linspace(t_u,t_c,10);
yy = feval(k_fun,tt);
plot(tt,yy)



% figure(1)
% plot(t,x(:,2),'r',t,x(:,3),'b',tm,Ibar(tm+1),'rx',tm,Rbar(tm+1),'bx');
% ax = gca;
% ax.XTick = t_0:14:t_f;
% ax.XTickLabel = date((t_0:14:t_f)+1);
% ax.XTickLabelRotation = 45;
% legend('I','R','I_{bar}','R_{bar}');
% title({
%     ['Dati Italia']
%     ['k_{c} = ' num2str(k_c,4)]
%     });
% 
% flag = x(:,1).*x(:,2)-2*beta*k_c;


% VINCOLO NON LINEARE PER MIN SU K

function [c,ceq] = mycon(K)

global x0 Nass tspan beta gamma
    
    SI = @(t,x) [-(beta - 0.5*x(1)*x(2)/K)*x(1)*x(2);
                  (beta - 0.5*x(1)*x(2)/K)*x(1)*x(2) - gamma*x(2)];
          
    Jac = @(t,x) [ -beta*x(2) + x(1)*(x(2)^2)/K, -beta*x(1) + (x(1)^2)*x(2)/K;
                    beta*x(2) - x(1)*(x(2)^2)/K,  beta*x(1) - (x(1)^2)*x(2)/K - gamma];
    options.Jacobian = Jac;
    
    [t, xm]  = eulerorosenbrock(SI,tspan,x0,options);
    xm(:,3) = ones(length(t),1) - xm(:,1) - xm(:,2);      % ricavo R per post-processing
    %xm = Nass.*xm;                                        % ???
    
    % controllo la condizione con i valori in percentuale, altrimenti e'
    % impossibile verificare la condizione
    
	% Nonlinear inequality constraints (c(K)<=0)
    c = xm(:,1).*xm(:,2) - 2*beta*K;   % deve essere <=0 (ATTENZIONE all'=)
    
    ceq = [];
end
