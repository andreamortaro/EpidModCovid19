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

global days K_disc regione ffit

if isempty(regione) % nel caso italia in questa funzione regione = []
    switch ffit
        case 0
            K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
        case 1
            K = @(t) A(1)*exp(-A(2)*t).*(1-exp(-A(3)*t)).^3;
    end
else
    %K = @(t) A(1)*exp(-A(2)*t).*(1-exp(-A(3)*t)).^3;
    K = @(t) A(1)*exp(-((t-A(2))/A(3)).^2);
end
    
% Calcolo numericamente la funzione dei minimi quadrati L
n = 2;
L = 0;
for jj = 1:1:length(days)
    day = days(jj);
    L = L+(K_disc(jj)-K(day)).^n;      % misura minimi quadrati
end

%deltat = days(end)-days(1);
%L = deltat*L;

%L = L/Nass;

end
