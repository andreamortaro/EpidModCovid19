%{
Prima del 09/03/20: il lockdown non e' ancora partito, l'obbiettivo ora e'
stimare i parametri beta e gamma del modello.
Aggiustamento versione 2020-06-01: analisi sensitivita'.
%}

%% Identificazione parametri:

close all
clear
clc

global  x0 tm ym Nass t_0 t_u pnt

% Upload dati protezione civile
tmp = fullfile('..','00 - dpc_data','2020-05-22','dati-andamento-nazionale');
[status,result] = fileattrib(tmp);
path_folder = result.Name;              % percorso alla cartella
[date,Ibar,Rbar] = data_read_dpc(path_folder);
% Attenzione: poiche' gli indici in matrice partono da 1 e tm parte da 0
% date(i), Ibar(i), Rbar(i) e' in corrispondenza con t_i-1

% DATI:
t_0 = 0;                            % 2020-02-24
t_u = 14;                           % 2020-03-09 t finale senza controllo
tm  = t_0:1:t_u;                    % tm = [0,..,t_u]
ym  = [Ibar(tm+1),Rbar(tm+1)];

K0  = [0.5,0.1];                    % guess iniziale per [beta,gamma]
pnt = 5;                            % piu nodi per migliore risoluzione sistema minquad

Nass = 60317000;                    % popolazione italiana istat 11.02.2020
I0 = Ibar(t_0+1); R0 = Rbar(t_0+1); S0 = Nass-I0-R0;
x0 = [S0;I0]/Nass;                  % dato iniziale in percentuale

% Minimizzazione

problem.options = optimoptions('fmincon','Display','iter');
problem.solver  = 'fmincon';
problem.objective = @minquad;
problem.x0 = K0;
problem.lb = [0,0];

K = fmincon(problem)    %#ok<NOPTS>
k1 = K(1); k2 = K(2);   % update parametri stimati 

SI = @(t,x) [-k1*x(1)*x(2); k1*x(1)*x(2) - k2*x(2)];   % update sys, nuovi parametri stimati

options.Jacobian = @(t,x) [-k1*x(2), -k1*x(2);
                            k1*x(2), k1*x(1) - k2];


[t, x]  = rk4(SI,[t_0,t_u],x0,options);     % sistema lineare

x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R per post-processing
x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass


% FIGURA

% imposto latex come inteprete per i grafici
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex');

fig = figure()
plot(t,x(:,2),'SeriesIndex',1);
hold on
plot(t,x(:,3),'SeriesIndex',2);
plot(tm,Ibar(tm+1),'o',...
    'MarkerSize',3,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[1 .6 .6]);
plot(tm,Rbar(tm+1),'o',...
    'MarkerSize',3,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);


ax = gca;
ax.XTick = 0:7:14;
ax.XTickLabel = date((0:7:14)+1);
ax.XTickLabelRotation = 45;
box on
legend('I','R','$I_{bar}$','$R_{bar}$','Location','NorthWest');
% title({
%     ['Dati Italia']
%     ['$\beta$ = ' num2str(k1,4)]
%     ['$\gamma$ = ' num2str(k2,3)]
%     });
ylabel('casi confermati');
set(gca,'FontSize',12.5)


exportgraphics(fig,'italia-preLock.pdf','ContentType','vector',...
               'BackgroundColor','none')

beta  = k1;
gamma = k2;
R_0   = beta/gamma;

table(beta,gamma,R_0)


%% Analisi sensitivita'

x1 = x(:,1)/Nass;	% S(t)
x2 = x(:,2)/Nass;  	% I(t)
x3 = x(:,3)/Nass;  	% R(t)

D0 = [0 0 0 0 0 0]';   % condizione iniziale del sistema D

% Matrice sensitivita'
% m non e' il tempo, ma l'indice che gira sul vettore x

sensys = @(m,D) [-beta*x2(m)*D(1)-beta*x1(m)*D(3)-x1(m)*x2(m);...
                 -beta*x2(m)*D(2)-beta*x1(m)*D(4);...
                  beta*x2(m)*D(1) + (beta*x1(m)-gamma)*D(3) + x1(m)*x2(m);...
                  beta*x2(m)*D(2) + (beta*x1(m)-gamma)*D(4) - x2(m);...
                  gamma*D(3);...
                  gamma*D(4) + x2(m)];
              
% Risolvo con Eulero Esplicito, il sistemare e' lineare in D
D = odesol(sensys,t,D0);

fig2 = figure();
plot(t,D(:,1),t,D(:,2),t,D(:,3),t,D(:,4),t,D(:,5),t,D(:,6))
%title('Sensitivity')
H=legend('$\partial_{\beta}S$','$\partial_{\gamma}S$',...
         '$\partial_{\beta}I$','$\partial_{\gamma}I$',...
         '$\partial_{\beta}R$','$\partial_{\gamma}R$');
set(H,'Location','Best',...
      'NumColumns',2,...
      'Orientation','horizontal')
axis tight
grid on
xlabel('t')
set(gca,'FontSize',12.5)

exportgraphics(fig2,'sensitivita.pdf','ContentType','vector',...
               'BackgroundColor','none')