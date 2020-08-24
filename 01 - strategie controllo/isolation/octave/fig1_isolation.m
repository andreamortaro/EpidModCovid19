% Modello SIR

clc
close all
clear all

% discretizzazione temporale:
T       = 20;                       % tempo finale
ntot    = 500;                      % numero di nodi
dt      = T/(ntot-1);               % passo
tspan   = linspace(0,T,ntot);

% parametri
beta    = 0.008;
mu      = 0.08;

N   = 100;
I0  = 20; S0 = N-I0; R0 = 0;

alfa = 1;
u_max = 80;

x0=[S0;I0];          % dato iniziale
options.InitialStep = dt;

% 1. CONTROLLO SU [0,T]

%sovrascrivo i t1,t2
t1=0;                % inizio controllo
t2=T;                % fine controllo

SIR = @(t,x) [-beta*x(1)*x(2);
               beta*x(1)*x(2)-mu*x(2)-alfa*x(2)*(t>=t1 & t<=t2)];

[tsol1,xsol1]= rk4(SIR,[0,T],x0,options);

% disegno dello scalino
scalino1 = zeros(ntot,1);
scalino1(tsol1>=t1 & tsol1<=t2) = u_max;
[xs1,ys1] = stairs(tsol1,scalino1);


% 2. CONTROLLO SU [before peak,T]

%sovrascrivo i t1,t2
t1=1;                % inizio controllo
t2=T;                % fine controllo

SIR = @(t,x) [-beta*x(1)*x(2);
               beta*x(1)*x(2)-mu*x(2)-alfa*x(2)*(t>=t1 & t<=t2)];

[tsol2,xsol2]= rk4(SIR,[0,T],x0,options);

% disegno dello scalino
scalino2 = zeros(ntot,1);
scalino2(tsol2>=t1 & tsol2<=t2) = u_max;
[xs2,ys2] = stairs(tsol2,scalino2);


% 3. CONTROLLO SU [time_peak,T]

%sovrascrivo i t1,t2
t1=4;                % inizio controllo
t2=T;                % fine controllo

SIR = @(t,x) [-beta*x(1)*x(2);
               beta*x(1)*x(2)-mu*x(2)-alfa*x(2)*(t>=t1 & t<=t2)];

[tsol3,xsol3]= rk4(SIR,[0,T],x0,options);

% disegno dello scalino
scalino3 = zeros(ntot,1);
scalino3(tsol3>=t1 & tsol3<=t2) = u_max;
[xs3,ys3] = stairs(tsol3,scalino3);


% 4. CONTROLLO SU [after peak,T]

%sovrascrivo i t1,t2
t1=10;               % inizio controllo
t2=T;                % fine controllo

SIR = @(t,x) [-beta*x(1)*x(2);
               beta*x(1)*x(2)-mu*x(2)-alfa*x(2)*(t>=t1 & t<=t2)];

[tsol4,xsol4]= rk4(SIR,[0,T],x0,options);

% disegno dello scalino
scalino4 = zeros(ntot,1);
scalino4(tsol4>=t1 & tsol4<=t2) = u_max;
[xs4,ys4] = stairs(tsol4,scalino4);


%% Stampo le figure

graphics_toolkit gnuplot;
colormap ("default");
set(gca, 'fontsize', 15);

spa = figure();
hax = plotyy (tsol1, xsol1(:,2),xs1,ys1);
set(hax,'ycolor','k','ylim',[0,N])
set(hax(2), 'YTick',[0,u_max])
set(hax(2), 'YTickLabel',{'0','$u^{max}$'})
axis([0 T 0 N]);
xlabel("t");
ylabel("I(t)");
grid on

spb = figure();
hax = plotyy (tsol2, xsol2(:,2),xs2,ys2);
set(hax,'ycolor','k','ylim',[0,N])
set(hax(2), 'YTick',[0,u_max])
set(hax(2), 'YTickLabel',{'0','$u^{max}$'})
xlabel("t");
ylabel("I(t)");
axis([0 T 0 N]);
grid on

spc = figure();
hax = plotyy (tsol3, xsol3(:,2),xs3,ys3);
set(hax,'ycolor','k','ylim',[0,N])
set(hax(2), 'YTick',[0,u_max])
set(hax(2), 'YTickLabel',{'0','$u^{max}$'})
xlabel("t");
ylabel("I(t)");
axis([0 T 0 N]);
grid on

spd = figure();
hax = plotyy (tsol4, xsol4(:,2),xs4,ys4);
set(hax,'ycolor','k','ylim',[0,N])
set(hax(2), 'YTick',[0,u_max])
set(hax(2), 'YTickLabel',{'0','$u^{max}$'})
xlabel("t");
ylabel("I(t)");
axis([0 T 0 N]);
grid on

% print dei vari subplot
print (spa, "subplot1", "-dpdflatexstandalone", "-F:20","-S640,480");
system ("pdflatex subplot1");
print (spb, "subplot2", "-dpdflatexstandalone", "-F:20","-S640,480");
system ("pdflatex subplot2");
print (spc, "subplot3", "-dpdflatexstandalone", "-F:20","-S640,480");
system ("pdflatex subplot3");
print (spd, "subplot4", "-dpdflatexstandalone", "-F:20","-S640,480");
system ("pdflatex subplot4");


