function vectorfield(deqns,xval,yval,t,N) 

%% Vector Field.
%
%   deqns   = sistema differenziale
%   xval    = tempi in cui valuto x
%   yval    = tempi in cui valuto y
%   N       = normalizzazione plot (se risolvo per [0,1]x[0,1] posso poi espandere a [0,N]x[0,N])

switch nargin
    case 3
        t=0;
    case 4
        N=1;
end

m=length(xval);
n=length(yval);
x1=zeros(n,m);
y1=zeros(n,m);

for a=1:m       % run colonne
    for b=1:n   % run righe
        pts = feval(deqns,t,[xval(a);yval(b)]);     % valuto il sistema diff in (t,[...])
        x1(b,a) = pts(1);
        y1(b,a) = pts(2);
    end
end

arrow = sqrt(x1.^2+y1.^2);      % modulo
quiver(xval*N,yval*N,x1./arrow,y1./arrow,.5,'r'); 
axis tight;
