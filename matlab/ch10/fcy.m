function fv = fcy(x)
% Corrugated spring function
C = 0; for j = 1:2, C = C + (x(j) - 5)^2; end
fv = -cos(5*sqrt(C)) + 0.1*C;
end
