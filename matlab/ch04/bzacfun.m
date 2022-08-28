function dy = bzacfun(x,y) 
% dy/dx for benzene/acetic acid system 
dy = (y*(1-y)/(y-x)) * (-1463.0572*x^3 + 2745.4788*x^2 - 1788.1398*x ...    
+ 556.8470)/(-365.7643*x^4 + 915.1596*x^3 - 894.0699*x^2 + 556.8470*x + 56.2570); 
end 