% vapcm.m
% differential equations are defined by the subfunction vapr
vapdat; % retrieve data
mL0 = 2800; mV0 = 100; z0 = [mV0 mL0]; tspan = [0 0.1]; [t z] = ode45(@vapr,tspan,z0);
plot(t,z(:,1),t,z(:,2),':'), grid, xlabel('t(h)'), ylabel('m_V,m_L(kg)'),
legend('m_V','m_L')


