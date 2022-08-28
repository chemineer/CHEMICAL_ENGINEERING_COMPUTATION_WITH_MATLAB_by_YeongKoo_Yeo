function [feedn, totaln] = bindistMT(alpha,q,zf,xd,xb,R)
% Calculation of binary distillation using McCabe/Thiele method
% input
%  alpha: relative volatility
%  q: feed condition parameter
%  zf,xd,xb: mole fractions of feed, distillate and bottoms
%  R: reflux ratio
% output
%  feedn: feed stage
%  totaln: number of total stages
%  Initialization and calculation of equilibrium curve
y = 0:0.1:1; ye = y; xe = ye./(alpha + (1-alpha)*ye);
xq = ((R+1)*zf + (q-1)*xd)/(R + q); yq = (R*zf + q*xd)/(R + q);
plot(xe,ye,'r'); hold on, axis([0 1 0 1]); set(line([0 1],[0 1]),'Color',[0 0 0]);
set(line([xd xq],[xd yq]),'Color',[1 0 1]); set(line([zf xq],[zf yq]),'Color',[1 0 1]);
set(line([xb xq],[xb yq]),'Color',[1 0 1]);
%Rectifying section
i = 1; xop(1) = xd; yop(1) = xd; y = xd;
while (xop(i) > xq)
xop(i+1) = y./(alpha + (1-alpha)*y);
yop(i+1)=R*xop(i+1)/(R+1) + xd/(R+1); % rectifying operating line
y = yop(i+1);
set(line([xop(i) xop(i+1)],[yop(i) yop(i)]),'Color',[0 0 1]);
if (xop(i+1) > xq), set(line([xop(i+1) xop(i+1)],[yop(i) yop(i+1)]),'Color', [0 0 1]); end
i = i+1;
end
feedn = i-1;
% Stripping section
c1 = (yq - xb)/(xq - xb); c2 = (yq - xq)/(xq - xb); yop(i) = c1*xop(i) - c2*xb; y = yop(i);
set(line([xop(i) xop(i)],[yop(i-1) yop(i)]),'Color',[0 0 1]);
while (xop(i)>xb)
xop(i+1) = y/(alpha + (1-alpha)*y); yop(i+1) = c1*xop(i+1) - c2*xb; y = yop(i+1);
set(line([xop(i) xop(i+1)],[yop(i) yop(i)]),'Color',[0 0 1]);
if (xop(i+1) > xb), set(line([xop(i+1) xop(i+1)],[yop(i) yop(i+1)]),'Color',[0 0 1]); end
i=i+1;
end
set(line([xop(i) xop(i)],[yop(i-1) yop(i)]),'Color',[0 0 1]); hold off,
xlabel('x'), ylabel('y'), totaln = i-1;
fprintf('Feed stage = %g\n', feedn); fprintf('Number of stages = %g\n', totaln);
end
