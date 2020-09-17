function pdecaller()
    options = odeset('OutputFcn',@odephas2);
    for i=-4:4
        for j= -4:4
              [t,y] = ode45(@vdp1,[0 20],[i;j],options); hold on;
        end
    end
end
function dydt = vdp1(t,y)
     dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];
end