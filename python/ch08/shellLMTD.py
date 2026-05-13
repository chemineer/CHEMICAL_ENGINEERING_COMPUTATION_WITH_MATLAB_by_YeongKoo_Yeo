import math

def shellLMTD():
    print("shellLMTD.m : number of required shells and LMTD")
    #T1 = float(input('Hot fluid inlet temperature (deg.F): '))
    #T2 = float(input('Hot fluid outlet temperature (deg.F): '))
    #t1 = float(input('Cold fluid inlet temperature (deg.F): '))
    #t2 = float(input('Cold fluid outlet temperature (deg.F): '))
    T1 = 250
    T2 = 100
    t1 = 80
    t2 = 120
    
    N, Dt1, Dt2 = 1, T2 - t1, T1 - t2
    P = (t2 - t1) / (T1 - t1)
    R = (T1 - T2) / (t2 - t1)
    
    A, B, F = -1.0, -1.0, 0.1
    
    while F <= 0.75:
        if R != 1:
            while A < 0:
                Pp = (1 - ((P * R - 1) / (P - 1))**(1/N)) / (R - ((P * R - 1) / (P - 1))**(1/N))
                A = (2 / Pp - 1 - R + math.sqrt(R**2 + 1)) / (2 / Pp - 1 - R - math.sqrt(R**2 + 1))
                if A < 0: N += 1
            F = (math.sqrt(R**2 + 1) * math.log10((1 - Pp) / (1 - Pp * R)) / (R - 1)) / math.log10(A)
        else: # R == 1
            while B < 0:
                Ppp = P / (N - P * (N - 1))
                B = (2 / Ppp - 2 + math.sqrt(2)) / (2 / Ppp - 2 - math.sqrt(2))
                if B < 0: N += 1
            F = (math.sqrt(R**2 + 1) * Ppp / (math.log(10 * (1 - Ppp)))) / math.log10(B)
            
        if F <= 0.75: N += 1
        
    LMTD = Dt1 if abs(Dt1 - Dt2) < 1e-9 else (Dt1 - Dt2) / math.log(Dt1 / Dt2)
    cLMTD = F * LMTD
    
    print(f'\nNumber of shells = {int(N)}')
    print(f'F factor = {F:9.4f}')
    print(f'Corrected LMTD = {cLMTD:9.4f}')