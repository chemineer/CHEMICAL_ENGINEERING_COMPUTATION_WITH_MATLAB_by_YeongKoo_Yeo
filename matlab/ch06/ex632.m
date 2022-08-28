% pmmarxn.m: PMMA polymerization reaction
% Data structure
pm.M0 = 1.5e4; pm.MWm = 0.10013; pm.MWi = 0.077; pm.Mjp = 0.18781;
pm.rhop = 1200; pm.Vms = 8.22e-4; pm.Vps = 7.7e-4; pm.Vis = 8.25e-4; pm.T = 350;
tspan = [0 6000]; M0 = pm.M0; x0 = zeros(1,10); x0(2) = M0; % Initial values
[t x] = ode15s(@pmmar,tspan,x0,[],pm);
I = x(:,1); M = x(:,2); R = x(:,3); L0 = x(:,4); L1 = x(:,5); L2 = x(:,6);
N0 = x(:,7); N1 = x(:,8); N2 = x(:,9); Q = x(:,10);
X = (M0 - M)/M0; % conversion
mom0 = L0 + N0; mom1 = L1 + N1; mom2 = L2 + N2;
smom0 = size(mom0); n = smom0(1,1);
for k = 2:n % calculate molecular weight
Mn(k,1) = mom1(k-1,1)/mom0(k-1,1); Mw(k,1) = mom2(k-1,1)/mom1(k-1,1);
end
for k = 1:n,  Pd(k,1) = Mw(k,1)/Mn(k,1); end
PMn = 100.13*Mn; PMw = 100.13*Mw;
% Plot results
subplot(3,2,1), plot(t,X), xlabel('t(s)'), ylabel('Conversion(X)')
subplot(3,2,2), plot(t,I), xlabel('t(s)'), ylabel('Initiator(moles)')
subplot(3,2,3), plot(t,M), xlabel('t(s)'), ylabel('Monomer(moles)')
subplot(3,2,4), plot(t,Q), xlabel('t(s)'), ylabel('Q(kJ)')
subplot(3,2,5), plot(t,PMn), xlabel('t(s)'), ylabel('Molecular weight PMn')
subplot(3,2,6), plot(t,PMw), xlabel('t(s)'), ylabel('Molecular weight PMw')
