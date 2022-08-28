% rxnt1dintp.m: performs 1D interpolation and plots the pchip interpolation 
tm = [0 18 42 55 70 82 86 95 102 115]; TC = [15 23 24 36 78 80 98 96 127 126]; % data 

tinv = linspace(0,115); yp = interp1(tm,TC,tinv,'pchip'); % 1D interpolation using the pchip method 
plot(tm,TC,'o',tinv,yp), xlabel('Time(min)'), ylabel('Temperature(deg.C)') 
axis([0 115 0 130]), legend('Data','Cubic Hermite interpolation','location','best')