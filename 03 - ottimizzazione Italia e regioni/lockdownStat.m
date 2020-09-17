function [tL, ImedioL, VarmediaL] = lockdownStat(data,I0f,R0f,beta_hist,gamma_hist,options)

% recupero i valori che mi servono
[Nass,~,~] = data.value;
[~,~,t_u,t_c,~] = data.time;
Kfun = data(3).Kvalue;

% settaggio per parametri
[alfab,alfag,p,M,B] = data.parametersStat;

deltatc = 0;
nstep = 500;
if isfield(options,'deltatc')
    deltatc = options.deltatc;
    nstep = options.nstep;
end

% salvo simulazioni (con rk4 non è facile prevedere la dimensione)
I_hist = cell(B,M); % per ogni riga beta,gamma fissato
t_hist = cell(B,M);
I_mean = cell(1,M); % valori medi, in ogni colonna ho un beta diverso


for ii = 1:B
    
    betaz = beta_hist(ii);
    gammaz = gamma_hist(ii);
    
    % update sistema
    SI = @(t,x) [-(betaz - x(1)*x(2)/Kfun(t))*x(1)*x(2);
                  (betaz - x(1)*x(2)/Kfun(t))*x(1)*x(2) - gammaz*x(2)];

    Jac = @(t,x) [ -betaz*x(2) + 2*x(1)*(x(2)^2)/Kfun(t), -betaz*x(1) + 2*(x(1)^2)*x(2)/Kfun(t);
                    betaz*x(2) - 2*x(1)*(x(2)^2)/Kfun(t),  betaz*x(1) - 2*(x(1)^2)*x(2)/Kfun(t) - gammaz];
    opt.Jacobian = Jac;
    
    %fig = figure('Visible','Off');
    fig = figure();
    hold on
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
        
    for jj = 1:M    % fissato il modello, simulo per diversi dati iniziali
        
        I0 = I0f(ii,jj); R0 = R0f(ii,jj); S0 = Nass - I0 - R0; 
        x0 = [S0;I0]/Nass;

        [t, x]  = eulerorosenbrock(SI,linspace(t_u,t_c+deltatc,nstep),x0,opt);
        x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);      % ricavo R per post-processing
        x = Nass.*x;                                       % normalizzo

        plot(t,x(:,2))
        legend('infetti')
        grid on
        box on
        
        % salvo la simulazione
        I_hist{ii,jj} = x(:,2);   % lo salvo in colonna
        t_hist{ii,jj} = t;
        
    end
    
    % fissato beta, calcolo la traiettoria I(t) media
    tmp = cell2mat(I_hist(ii,:));
    mean = sum(tmp,2)./M;
    I_mean{1,ii} = mean;

    hold on
    A = plot(t,mean,'LineStyle','--','LineWidth',1.5,'Color','black');
    legend(A,'media')

    txt = sprintf('I0=%d, M =%i, p=%.2f',I0,M,p);
    title(txt)

    set(gca,'FontSize',12.5);
end

% media empirica
tmp = cell2mat(I_mean);
Imedio = sum(tmp,2)./B;
%h = plot(t,Imedio,'black','LineWidth',3,'Linestyle','--');
%legend(h,'$\mu$')


%% Calcolo Varianza
% per ogni tempo tk calcolo la varianza, cioe' distanza media dei valori I(tk)
% dal valore medio Imedio

% Welford
Var = cell(1,B);
% per ogni beta fissato calcolo la varianza delle varie curve I(t)
for ii = 1:B
    % traiettorie per beta fissato
    tmp = cell2mat(I_hist(ii,:));
    tmpVar = zeros(length(t),1);
    for tk = 1:1:length(t)
        %tmpVar(tk) = var(tmp(tk,:));
        tmpVar(tk) = online_variance(tmp(tk,:));
    end
    Var{1,ii} = tmpVar;
    
    figure('Visible','off');
    %figure();
    plot(t_hist{ii,1},tmpVar)
end
% calcolo una varianza media tra le varianze ottenute per diversi beta
tmp = cell2mat(Var);
Varmedia = sum(tmp,2)./B;

% confronto tra calcolo e funzione built-in
figure();
plot(t,Varmedia)
xlabel('t')
title('Varianza media')
legend('welford')
grid on
set(gca,'FontSize',12.5);


% salvo i dati
tL = t_hist{1,1};
ImedioL = Imedio;
VarmediaL = Varmedia;

end