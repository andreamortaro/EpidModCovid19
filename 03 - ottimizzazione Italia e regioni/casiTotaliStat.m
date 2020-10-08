function [cTmedio,cTVarmedia,cT] = casiTotaliStat(data,tt,hist,ssave)

%
%   [cTmedio,cTVarmedia,cT] = casiTotaliStat(data,tt,hist,ssave)
%
%   Date le simulazioni degli infetti e rimossi nella struttura hist, la
%   funzione casiTotaliStat calcola i casi totali cT e restituisce i casi 
%   totali medi cTmedio, la loro varianza cTVarmedia.
%

% recupero i dati dalla struttura data
[t_0,~,t_u,t_c,date] = data.time;
[~,~,~,totCases] = data.value;
[~,~,~,M,B] = data.parametersStat;

% recupero le simulazioni
I_hist = hist(5).sim; R_hist = hist(6).sim;
casiTotali = cell(B,M);     % salvo i casi totali per ogni simulazione
cT_mean = cell(1,B);        % salvo i casi totali medi per ogni (beta,gamma) fissato
for ii = 1:B
    for jj = 1:M
        II = I_hist{ii,jj}; RR = R_hist{ii,jj};     % fisso una simulazione
        casiTotali{ii,jj} = RR+II;
    end
    % fissato beta, calcolo la traiettoria dei casi totali media
    tmp = cell2mat(casiTotali(ii,:));
    mean = sum(tmp,2)./M;
    cT_mean{1,ii} = mean;
end

% media empirica
tmp = cell2mat(cT_mean);
cTmedio = sum(tmp,2)./B;
len = length(cTmedio);

% Welford
Var = cell(1,B);
% per ogni beta fissato calcolo la varianza delle varie curve I(t)
for ii = 1:B
    % traiettorie per beta fissato
    tmp = cell2mat(casiTotali(ii,:));
    tmpVar = zeros(len,1);
    for tk = 1:1:len
        %tmpVar(tk) = var(tmp(tk,:));
        tmpVar(tk) = online_variance(tmp(tk,:));
    end
    Var{1,ii} = tmpVar;
    
    figure('Visible','off');
    %figure();
    plot(tt,tmpVar)
end
% calcolo una varianza media tra le varianze ottenute per diversi beta
tmp = cell2mat(Var);
cTVarmedia = sum(tmp,2)./B;

%% figura

cT = cell2mat(casiTotali(:)');

ctfig = figure();

plot(t_0:t_c,totCases,'o',...
    'MarkerSize',4,...
    'MarkerEdgeColor',[.9 .5 .1],...
    'MarkerFaceColor',[.9 .7 .3]);
      
plot(tt,cTmedio,'Linewidth',2.5,'color','black'); p2.Color(4) = 0.6;

% grafico con intervallo di confidenza

% calcolo l'intervallo di confidenza al 95%
p = 0.05; CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CI95cT = zeros(length(tt),2);
for kk = 1:length(tt)
    x = cT(kk,:);
    CII = CIFcn(x,100-p*100);
    CI95cT(kk,:) = CII;
end

% calcolo l'intervallo di confidenza al 50%
p = 0.50; CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CI50cT = zeros(length(tt),2);
for kk = 1:length(tt)
    x = cT(kk,:);
    CII = CIFcn(x,100-p*100);
    CI50cT(kk,:) = CII;
end

plot(t_0:t_c,totCases,'o',...
    'MarkerSize',4,...
    'MarkerEdgeColor',[.9 .5 .1],...
    'MarkerFaceColor',[.9 .7 .3]);

hold on
H2 = plot(tt, CI95cT(:,1),...
          tt, CI95cT(:,2),'Color', [.5 .8 .3],'Visible','Off');
H3 = plot(tt, CI50cT(:,1),...
          tt, CI50cT(:,2),'Color', [.5 .7 .2],'Visible','Off');
ttt = [tt', fliplr(tt')]';
yyy = [H2(1).YData, fliplr(H2(2).YData)]'; fill(ttt,yyy,[.5 .8 .3],'facealpha',0.5,'LineStyle','none');
yyy = [H3(1).YData, fliplr(H3(2).YData)]'; fill(ttt,yyy,[.5 .7 .2],'facealpha',0.6,'LineStyle','none');

H1 = plot(tt, cTmedio, 'Color', [.3 .5 0], 'LineWidth', 2);

ax = gca;
ax.XTick = [t_u,37,67,98,t_c];
ax.XTickLabel = date([t_u,37,67,98,t_c]+1);
ax.XTickLabelRotation = 45;
xline(t_u,':','inizio Lockdown')
ylabel('casi totali attesi');

grid off
box on

if exist('regione','var') == 1
    title(char(regione));
else
    title('Italia')
end

set(gca,'FontSize',12.5);

%ctfig.CurrentAxes.XLim = [0 157];

if ssave == 1
    if exist('regione','var') == 1
        exportgraphics(ctfig,...
        ['figure/' char(regione) '/CTriepilogoStat ' num2str(date(t_c+1)) '.pdf'],...
        'ContentType','vector',...
        'BackgroundColor','none')
    else
        exportgraphics(ctfig,['figure/Italia/CTriepilogoStat ' num2str(date(t_c+1)) '.pdf'],...
         'ContentType','vector',...
        'BackgroundColor','none')
    end
end

end