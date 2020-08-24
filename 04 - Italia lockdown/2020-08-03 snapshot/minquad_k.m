function L = minquad_k(A)     % funzionale da minimizzare

%
%   L = minquad(K)
%
%   Definisco il funzionale L da minimizzare
%
%   INPUTS:
%   K       : incognita
%
%   OUTPUTS:
%   L       : funzione minimi quadrati
%

global  k_history

k_fun = @(t) A(1)*exp(A(2)*t)*(1-exp(A(3)*t));

% Calcolo numericamente la funzione dei minimi quadrati L
%psi = 0.5;
%phi = 1-psi;
n = 2;
iter = k_history(:,1);
days = k_history(:,2);
k_c  = k_history(:,3);         % recupero tutti i k stimati

L = 0;
for j = iter'               % traspongo per poter iterare sulle colonne (matlab)
    day = days(j);
    L = L+(k_c(j)-k_fun(day)).^n;
%     L = L+(psi*(ym(j,1)-xm(tt(j),2)).^n + phi*(ym(j,2)-xm(tt(j),3)).^n);  % misura minimi quadrati
end

deltat = days(end)-days(1);
L = deltat*L;

end
