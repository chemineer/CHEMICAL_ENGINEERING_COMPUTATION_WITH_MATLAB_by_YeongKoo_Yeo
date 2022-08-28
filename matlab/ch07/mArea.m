function y = mArea(x,D,E,F,R,S,T,alpa,r,xf,uf)
u = -D.*x + sqrt((D.^2).*(x.^2) + 2*E.*x + F.^2);
fi = (D.*x - F) + sqrt((D.^2)*(x.^2) + 2*E.*x + F.^2);
ud = (1-xf).*((uf-E/D)./(u-E/D)).^R.*((uf-alpa+F)./(u-alpa+F)).^S...
.*((uf-F)./(u-F)).^T;
y = ud./((fi-x).*(1./(1+x) - r./(1+fi)));
end
