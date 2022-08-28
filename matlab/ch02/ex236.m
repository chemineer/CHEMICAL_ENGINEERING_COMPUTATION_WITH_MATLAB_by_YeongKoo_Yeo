% hum2dintp.m: performs 2D interpolation 
RH  =  [10 30 50 70 90]; T  =  [51 44]; 
H  =  [8.27 24.75 41.30 58.10 73.51; 6.52 19.58 32.70 45.75 58.78]; 
DP  =  [10.10 26.89 37.01 43.21 47.80; 6.40 23.34 32.18 38.32 42.81]; 
Hv  =  interp2(RH,T,H,58.4,46.8,'spline') % absolute humidity 
DPv  =  interp2(RH,T,DP,58.4,46.8,'spline') % dew point