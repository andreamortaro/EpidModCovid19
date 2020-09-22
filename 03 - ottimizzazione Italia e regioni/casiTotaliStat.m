function [cTmedio,cTVarmedia] = casiTotaliStat(data,tt,hist)

[~,~,~,M,B] = data.parametersStat;

[I_hist,R_hist] = hist.sim;
casiTotali = cell(B,M);
cT_mean = cell(1,B);
for ii = 1:B
    for jj = 1:M
        II = I_hist{ii,jj}; RR = R_hist{ii,jj};
        casiTotali{ii,jj} = II + RR;
    end
    % fissato beta, calcolo la traiettoria I(t) media
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

end

