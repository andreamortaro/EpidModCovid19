function L = minquad_kcontinuo(A)     % funzionale da minimizzare

%
%   L = minquad_kcontinuo(A)
%
%   Definisco il funzionale L da minimizzare.
%
%   INPUTS:
%   A       : incognita della forma A = [A(1),A(2),A(3)]
%
%   OUTPUTS:
%   L       : funzione minimi quadrati
%

global days K_disc t_c

switch t_c
    case 127
        K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
    case 188
        K = @(t) A(1)*exp(-A(2)*t).*(1-exp(-A(3)*t)).^2;
%         K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
    otherwise
        K = @(t) -A(1)*t^2 + A(2)*t - A(3);
end


% Calcolo numericamente la funzione dei minimi quadrati L
n = 2;
L = 0;
for jj = 1:1:length(days)
    day = days(jj);
    L = L+(K_disc(jj)-K(day)).^n;      % misura minimi quadrati
end

deltat = days(end)-days(1);
L = deltat*L;

end