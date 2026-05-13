import numpy as np

def brentopt(objfun, x1, h, crit):
    """
    Brent's algorithm for 1-dimensional minimization.
    """
    clarge = 1e40
    tau = (np.sqrt(5) - 1) / 2
    nf = 1
    f1 = objfun(x1)
    
    # Determine initial 3-points
    x2 = x1 + h
    nf += 1
    f2 = objfun(x2)
    
    if f2 > f1:
        x1, x2 = x2, x1
        f1, f2 = f2, f1
        h = -h
        
    # Search phase
    while x1 < clarge:
        h = h / tau
        x4 = x2 + h
        nf += 1
        f4 = objfun(x4)
        if f4 > f2:
            break
        f1, x1 = f2, x2
        f2, x2 = f4, x4
        
    x3 = x4
    f3 = f4
    a, b = (x1, x3) if x1 < x3 else (x3, x1)
    fa, fb = (f1, f3) if x1 < x3 else (f3, f1)
    
    x = x2
    fx = f2
    w = v = x
    fw = fv = fx
    ev = 0
    d = 0 # 추가된 초기화
    
    while True:
        hm = (a + b) / 2
        # Check interval convergence
        if abs(x - hm) <= crit - (b - a) / 2:
            return x, fx, nf
        
        etemp = 0
        p = 0
        q = 0
        
        if abs(ev) > crit:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2 * (q - r)
            if q > 0: p = -p
            else: q = -q
            etemp = ev
            ev = d
            
        ind1 = -1 if (q * (a - x) - p) < 0 else 1
        ind2 = -1 if (q * (b - x) - p) < 0 else 1
        
        if (abs(p) >= abs(q * etemp / 2)) or (ind1 == ind2):
            ev = (b - x) if x < hm else (a - x)
            d = (1 - tau) * ev
        else:
            d = p / q
            u = x + d
            if (u - a < crit) or (b - u < crit):
                d = crit if x < hm else -crit
        
        if abs(d) >= crit:
            u = x + d
        else:
            u = x + crit if d > 0 else x - crit
            
        nf += 1
        fu = objfun(u)
        
        # Update logic
        if fu <= fx:
            if u < x: b = x; fb = fx
            else: a = x; fa = fx
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else:
            if u < x: a = u; fa = fu
            else: b = u; fb = fu
            
            if (fu <= fw) or (w == x):
                v, fv = w, fw
                w, fw = u, fu
            elif (fu <= fv) or (v == x) or (v == w):
                v, fv = u, fu
                
        # Check function convergence
        if (fa - fx) + (fb - fx) < crit:
            return x, fx, nf