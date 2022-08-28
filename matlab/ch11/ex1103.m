% compdfz.m
x = -5:0.1:5;
A = trimf(x,[-5 -4 -2]); B = trapmf(x,[-5 -3 2 5]); C = max(0.7*A, 0.5*B); %produce fuzzy sets
plot(x,C,'LineWidth',4), ylim([-1 1])
x1 = defuzz(x,C,'centroid')
h1 = line([x1 x1],[-0.2 1.2],'Color','k'); t1 = text(x1,-0.2,'centroid','FontWeight','bold');
x2 = defuzz(x,C,'bisector')
h2 = line([x2 x2],[-0.4 1.2],'Color','k'); t2 = text(x2,-0.4,'bisector','FontWeight','bold');
x3 = defuzz(x,C,'mom')
h3 = line([x3 x3],[-0.7 1.2],'Color','k'); t3 = text(x3,-0.7,'mom','FontWeight','bold');
