function [tPL, ImedioPL, VarmediaPL,I0f,R0f,hist] = prelockStat(data,options)

%
%   prelockStat simula il modello SIR M*B volte e stampa a video
%   i casi positivi attesi durante il periodo antecedente al lockdown.
%   Nell'output della funzione salva i risultati ottenuti e le simulazioni
%   effettuate.
%

% recupero i valori che mi servono
[Nass,Ibar,Rbar] = data.value;
[t_0,~,t_u,~,date] = data.time;
[beta,gamma] = data.parameters;

% settaggio per parametri
[pl,~] = data.infSim;       % simulazione curva infetti in prelock e lock
[alfab,alfag,p,M,B] = data.parametersStat;

I0m = Ibar(t_0+1); R0m = Rbar(t_0+1);   % valori misurati

% se non definisco options
if nargin == 1
    ffig = 1;
else
    if isfield(options,'ffig')
        ffig = options.ffig;
    end
end

if isfield(data,'regione')
    regione = data(1).regione;
end

% salvo simulazioni (con rk4 non è facile prevedere la dimensione)
I_hist = cell(B,M); % per ogni riga beta,gamma fissato
R_hist = cell(B,M);
t_hist = cell(B,M);

%per il lockdown
beta_hist = zeros(B,1);
gamma_hist = zeros(B,1);
I0f = zeros(B,M);           % salvo i valori finali
R0f = zeros(B,M);
I_mean = cell(1,B);         % valori medi, in ogni colonna ho un beta diverso

for ii = 1:B

    % definisco i parametri betaz e gammaz
    z2 = random('Beta',2,2);
    %z2 = rand;
    betaz = beta + alfab*z2;
    gammaz = gamma + alfag*z2;
    
    % definisco il modello
    SI = @(t,x) [-betaz*x(1)*x(2);...
                  betaz*x(1)*x(2) - gammaz*x(2)];

    Jac = @(t,x) [-betaz*x(2), -betaz*x(2);
                   betaz*x(2),  betaz*x(1) - gammaz];
    opt.Jacobian = Jac;
    
    figure('Visible','Off');
    %figure();
    hold on
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
        
    for jj = 1:M    % fissato il modello, simulo per diversi dati iniziali
        
        % dato iniziale simulazione
        z1 = random('Beta',10,10);
        %z1 = rand;
        
        I0 = I0m*(1+p*z1); R0 = R0m*(1+p*z1); S0 = Nass - I0 - R0; 
        x0 = [S0;I0]/Nass;
        
        [t, x]  = rk4(SI,[t_0,t_u],x0,opt); 
        x(:,3) = ones(length(t),1) - x(:,1) - x(:,2);   % ricavo R per post-processing
        x = Nass.*x;                                    % ri-normalizzo da percentuale a Nass

        plot(t,x(:,2))
        legend('infetti')
        grid on
        box on
        
        % salvo la simulazione
        I_hist{ii,jj} = x(:,2);   % lo salvo in colonna
        R_hist{ii,jj} = x(:,3);   % lo salvo in colonna
        t_hist{ii,jj} = t;
        
    end

    beta_hist(ii) = betaz;
    gamma_hist(ii) = gammaz;
    
    tmp = cell2mat(R_hist(ii,:));
    R0f(ii,:) = tmp(end,:);
    %R_mean{1,ii} = sum(tmp,2)./M;
    
    % fissato beta, calcolo la traiettoria I(t) media
    tmp = cell2mat(I_hist(ii,:));
    mean = sum(tmp,2)./M;
    I_mean{1,ii} = mean;
    I0f(ii,:) = tmp(end,:);

    hold on
    A = plot(t,mean,'LineStyle','--','LineWidth',1.5,'Color','black');
    legend(A,'media')
end

% media empirica
tmp = cell2mat(I_mean);
Imedio = sum(tmp,2)./B;

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
        tmpVar(tk) = var(tmp(tk,:));
        %tmpVar(tk) = online_variance(tmp(tk,:));
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
%figure();
figure('Visible','off');
plot(t,Varmedia)
xlabel('t')
title('Varianza media - prelock')
legend('welford')
grid on
set(gca,'FontSize',12.5);

%% plot media e varianza

if ffig == 1
    
    fig = figure();
    hold on;      % plot media e varianza

    % grafico infetti PC
    
    % imposto latex come inteprete per i grafici
    set(groot,...
        'defaulttextinterpreter','latex',...
        'defaultAxesTickLabelInterpreter','latex',...
        'defaultLegendInterpreter','latex');
    
    tm = t_0:1:t_u;
    
    p1 = plot(tm,Ibar(tm+1),'o',...
        'MarkerSize',4,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6]);        
    hold on
    p3 = plot(pl(:,1),pl(:,2),'SeriesIndex',2,'Linewidth',2.5); p3.Color(4) = 0.6;
    

    % grafico con intervallo di confidenza
    
    %sm = sqrt(Varmedia);                % scarto quadratico medio
    
    p = 0.05; CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
    CI95 = zeros(length(t),2);
    tmp = cell2mat(I_hist(:)');
    for kk = 1:length(t)
        x = tmp(kk,:);
        CII = CIFcn(x,100-p*100);
        CI95(kk,:) = CII;
    end

    p = 0.50; CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
    CI50 = zeros(length(t),2);
    tmp = cell2mat(I_hist(:)');
    for kk = 1:length(t)
        x = tmp(kk,:);
        CII = CIFcn(x,100-p*100);
        CI50(kk,:) = CII;
    end


    hold on
    H2 = plot(t, CI95(:,1),...
              t, CI95(:,2),'Color', [.4 .5 .8],'Visible','Off');
    H3 = plot(t, CI50(:,1),...
              t, CI50(:,2),'Color', [.6 .8 1],'Visible','Off');
    tt = [t', fliplr(t')]';
    yy = [H2(1).YData, fliplr(H2(2).YData)]'; fill(tt,yy,[.4 .5 .8],'facealpha',0.5);
    yy = [H3(1).YData, fliplr(H3(2).YData)]'; fill(tt,yy,[.6 .8 1],'facealpha',0.5);
    
    grid on
    box on
    H1 = plot(t, Imedio, 'Color', 'k', 'LineWidth', 2);
    legend([H1, H2(1),H3(1)], '$\mu$', '0.95\%','0.50\%','Location', 'Best');


    ax = gca;
    ax.XTick = 0:7:14;
    ax.XTickLabel = date((0:7:14)+1);
    ax.XTickLabelRotation = 45;
    box on
    ylabel('casi attuali attesi');
    
end

% salvo i dati
tPL = t_hist{1,1};
ImedioPL = Imedio;
VarmediaPL = Varmedia;

hist(1).parameters = beta_hist;
hist(2).parameters = gamma_hist;
hist(1).sim = I_hist;
hist(2).sim = R_hist;

end
