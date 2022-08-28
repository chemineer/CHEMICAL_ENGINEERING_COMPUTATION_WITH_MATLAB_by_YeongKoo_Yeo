[ppnet,tr] = train(ppnet,x,t); 
y = ppnet(x); % output calculation by network 
neterr = gsubtract(y,t); % error 
perm = perform(ppnet,t,y) % performance 