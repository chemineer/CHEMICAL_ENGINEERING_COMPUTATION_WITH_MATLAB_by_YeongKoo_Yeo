%Conjugate Gradients Method 
A = [4 -1 -1 0; -1 4 0 -1; -1 0 4 -1; 0 -1 -1 4]; b = [45 35 55 45]'; 
x = pcg(A,b)