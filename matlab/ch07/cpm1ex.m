function res = cpm1ex
cpmdata;
Pa = Pm(1); a = 1 - alpa; b = -1 + alpa + 1./r + (alpa-1)*xr./r;
c = -alpa*xr./r; yp = (-b+sqrt(b.^2 - 4*a.*c))./(2*a); % permeate mole fraction
theta = (xf-xr)./(yp-xr); % stage-cut
Am = (theta.*qf.*yp)./((Pa./t).*(ph.*xr - pl.*yp)); % membrane area
rc = qf*theta*yp/(qf*xf); % recovery ratio
xom = (xf*(1+r*(alpa-1)*(1-xf)))/(xf*(1-alpa)+alpa);
% Results: res = [yp, xr, Am, theta, rc]
res.yp = yp; res.theta = theta; res.Am = Am; res.rc = rc; res.xr = xr;
end
