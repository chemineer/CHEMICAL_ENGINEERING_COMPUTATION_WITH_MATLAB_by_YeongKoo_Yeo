%Cracking of Acetone in a Plug-Flow Reactor

pf = [162 0]; Vspan = [0 4]; X0 = [38.3 0 0 1150]; [V X] = ode45(@adfun, Vspan, X0, [], pf);
fprintf('\nFA = %g,  FB = %g,  FC = %g, T  = %g\n', X(end,1), X(end,2), X(end,3), X(end,4))
