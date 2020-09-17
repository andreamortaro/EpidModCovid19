function [tPL, ImedioPL, VarmediaPL,I0f,R0f,beta_hist,gamma_hist] = prelockStat(data,options)


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
    ssave = 1;
else
    if isfield(options,'ffig')
        ffig = options.ffig;
    end
    if isfield(options,'ssave')
        ssave = options.ssave;
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
I_mean = cell(1,B); % valori medi, in ogni colonna ho un beta diverso
%R_mean = cell(1,B);

for ii = 1:B

    % definisco i parametri betaz e gammaz
    %z2 = random('Beta',2,2);
    z2 = rand;
    % beta e gamma come nell'articolo, ma così gamma diventa enorme
    betaz = beta + alfab*z2;
    gammaz = gamma + alfag*z2;
    
    % definisco il modello
    SI = @(t,x) [-betaz*x(1)*x(2);...
                  betaz*x(1)*x(2) - gammaz*x(2)];

    Jac = @(t,x) [-betaz*x(2), -betaz*x(2);
                   betaz*x(2),  betaz*x(1) - gammaz];
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

    beta_hist(ii) = beta;
    gamma_hist(ii) = gamma;
    
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
    
    txt = sprintf('I0=%d, M =%i, p=%.2f prelock ',I0,M,p);
    
    if exist('regione','var') == 1
        title([txt char(regione)]);
    else
        title([txt 'Italia'])
        
    end
    
%     % controllare se cell2mat(t_hist(1,:))
%     tmp = cell2mat(t_hist(1,:));
%     for kk = 1:length(t)
%         if isequal(tmp(kk,:)) == 0
%             warning('Attenzione ai tempi della soluzione')
%         end
%     end

end

% media empirica
tmp = cell2mat(I_mean);
Imedio = sum(tmp,2)./B;

% tmp = cell2mat(R_mean);
% Rmedio = sum(tmp,2)./B;

% % come nell'articolo
% mm = I0 + I0*p*.5;
% mm/I0

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
    
    sm = sqrt(Varmedia);                % scarto quadratico medio

    z_alpha1 = 1.96;                    % alpha = 0.05 (confidenza) --> liv conf 0.95
    delta1 = z_alpha1*sm/sqrt(M*B);     % giusto M*B?

    z_alpha2 = 1.2816;                  % alpha = 0.20 (confidenza) --> liv conf 0.80
    delta2 = z_alpha2*sm/sqrt(M*B);

    z_alpha3 = 0.67;                    % liv confidenza 0.50
    delta3 = z_alpha3*sm/sqrt(M*B);

    H1 = plot(t, Imedio, 'Color', 'k', 'LineWidth', 2);
    hold on
    H2 = plot(t, Imedio - delta1,...
              t, Imedio + delta1,'Color', 'b');
    H3 = plot(t, Imedio - delta3,...
              t, Imedio + delta3,'Color', 'b');
    tt = [t', fliplr(t')]';
    hold on
    yy = [H2(1).YData, fliplr(H2(2).YData)]'; h2 = fill(tt,yy,'b','facealpha',0.2);
    yy = [H3(1).YData, fliplr(H3(2).YData)]'; h3 = fill(tt,yy,'cyan','facealpha',0.2);
    legend([H1, H2(1),H3(1)], '$\mu$', '0.95\%','0.50\%','Location', 'Best');
    grid on
    box on

    ax = gca;
    ax.XTick = 0:7:14;
    ax.XTickLabel = date((0:7:14)+1);
    ax.XTickLabelRotation = 45;
    box on
    %legend([p1,p2,p3,p4],'$I_{bar}$','$R_{bar}$','I','R','Location','NorthWest');
    ylabel('casi confermati');
    
    if exist('regione','var') == 1
        title(char(regione));
    else
        title('Italia')
    end
    
    set(gca,'FontSize',12.5)
    if ssave == 1
        if exist('regione','var') == 1
            exportgraphics(fig,'figure/' + regione + '/prelockStat.pdf',...
            'ContentType','vector',...
            'BackgroundColor','none')
        else
            exportgraphics(fig,'figure/Italia/Italia_prelockStat.pdf',...
            'ContentType','vector',...
            'BackgroundColor','none')
        end
    end
end

% salvo i dati
tPL = t_hist{1,1};
ImedioPL = Imedio;
VarmediaPL = Varmedia;


end