function ind = compID(cname) 
% Assignment of identification number 
cname = upper(cname);  
switch cname    
case {'FLUORINE','F2'}, ind = 1;                
case {'CHLORINE','CL2'}, ind = 2;    
case {'SULFUR DIOXIDE','SO2'}, ind = 3;       
case {'CARBON MONOXIDE','CO'}, ind = 4;  
case {'CARBON DIOXIDE','CO2'}, ind = 5;       
case {'HYDROGEN CHLORIDE','HCL'}, ind = 6;     
case {'AMMONIA','NH3'}, ind = 7;                   
case {'WATER','H2O'}, ind = 8;     
case {'HYDROGEN PEROXIDE','H2O2'}, ind = 9; 
case {'HYDROGEN','H2'}, ind = 10;    
case {'NITROGEN','N2'}, ind = 11;                 
case {'OXYGEN','O2'}, ind = 12;    
case {'ETHYLENE','C2H4'}, ind = 13;             
case {'METHANE','CH4'}, ind = 14;    
case {'ETHANE','C2H6'}, ind = 15;                
case {'PROPANE','C3H8'}, ind = 16;    
case {'BENZENE','C6H6'}, ind = 17;              
case {'TOLUENE','C7H8'}, ind = 18;    
case {'ANILINE','C6H7N'}, ind = 19;            
case {'PHENOL','C6H6O'}, ind = 20;    
case {'CYCLOPROPANE','C3H6'}, ind = 21;     
case {'CYCLOHEXANE','C6H12'}, ind = 22;     
case {'1,3 BUTADIENE','C4H6'}, ind = 23;    
case {'METHANOL','CH4O'}, ind = 24;    
case {'CHLOROFORM','CHCL3'}, ind = 25;    
case {'CARBON TETRACHLORIDE','CCL4'}, ind = 26; 
end 