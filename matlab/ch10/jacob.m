function Hs = jacob(fun, x, h)
% Jacobian of f(x)
hd = 2*h; n = length(x); x = x(:)'; M = eye(n);
for k = 1:n, Hs(:,k) = (fun(x + M(k,:)*h) - fun(x-M(k,:)*h))'/hd; end
end
