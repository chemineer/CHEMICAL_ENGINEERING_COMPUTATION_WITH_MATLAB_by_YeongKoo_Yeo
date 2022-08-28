% muffs.m
% Create an FIS object
rxnfis = mamfis('Name',"cwrate");
% Add the first input variable (color) and membership functions
rxnfis = addInput(rxnfis,[0 10],'Name',"color");
rxnfis = addMF(rxnfis,"color","gaussmf",[1.5 0],'Name',"red");
rxnfis = addMF(rxnfis,"color","gaussmf",[1.5 5],'Name',"orange");
rxnfis = addMF(rxnfis,"color","gaussmf",[1.5 10],'Name',"yellow");
% Add the second input variable (pressure) and trapezoidal membership functions
rxnfis = addInput(rxnfis,[0 10],'Name',"pressure");
rxnfis = addMF(rxnfis,"pressure","trapmf",[-2 0 1 3],'Name',"low");
rxnfis = addMF(rxnfis,"pressure","trapmf",[7 9 10 12],'Name',"high");
% Add the output variable (coolant) and triangular membership functions
rxnfis = addOutput(rxnfis,[0 40],'Name',"coolant");
rxnfis = addMF(rxnfis,"coolant","trapmf",[-2 0 5 10],'Name',"low");
rxnfis = addMF(rxnfis,"coolant","trimf",[10 15 20],'Name',"average");
rxnfis = addMF(rxnfis,"coolant","trapmf",[20 25 40 50],'Name',"high");
% Create and add a numeric array representing rules
numrule = [1 1 1 1 2; 2 0 2 1 1; 3 2 3 1 2];
rxnfis = addRule(rxnfis,numrule);
% Evaluation of the created FIS object
in1 = evalfis(rxnfis,[1 2])
vals = [3 5; 2 7; 3 1]; % three different types of input variables
in3 = evalfis(rxnfis,vals)
gensurf(rxnfis)
