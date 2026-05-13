import math

def calculate_htcoef():
    # --- Data (입력 데이터) ---
    # Tube-side (관측)
    Fft = 0.002
    Kt = 0.366
    mut = 0.72
    Cpt = 1
    Jt = 9.75
    
    # Shell-side (동측)
    Ffs = 0.002
    Ks = 0.074
    mus = 2.45
    Cps = 0.518
    
    # Dimensions (크기 - ft 단위 변환 포함)
    Dis = 3.068 / 12  # shell inside diameter
    Dos = 3.5 / 12    # shell outside diameter
    Dit = 1.61 / 12   # tube inside diameter
    Dot = 1.9 / 12    # tube outside diameter
    
    # Fin size (핀 규격)
    N = 24
    Hf = 0.5 / 12
    w = 0.035 / 12
    
    # Flow rates (유량)
    Gt = 26740
    Gs = 18000
    
    # --- Tube-side Calculations (관측 계산) ---
    Deq_t = Dit
    Pr_t = Cpt * (2.4191 * mut) / Kt  # 1 cP = 2.4191 lb/(ft-hr)
    Re_t = Gt * Deq_t / (2.4191 * mut)
    
    if Re_t < 10000:
        hit = Jt * (Pr_t**(1/3)) * Kt / Deq_t
    else:
        hit = 0.023 * (Re_t**0.8) * (Pr_t**(1/3)) * Kt / Deq_t
        
    hift = 1 / (1 / hit + Fft)
    
    # --- Shell-side (bare tube) Calculations (핀이 없는 상태의 동측 계산) ---
    Deq_s_bare = Dis - Dot
    Pr_s = Cps * (2.4191 * mus) / Ks
    Re_s = Gs * Deq_s_bare / (2.4191 * mus)
    
    if Re_s < 10000:
        his = Jt * (Pr_s**(1/3)) * Ks / Deq_s_bare
    else:
        his = 0.023 * (Re_s**0.8) * (Pr_s**(1/3)) * Ks / Deq_s_bare
        
    hifs = 1 / (1 / his + Ffs)
    
    # --- Shell-side (finned-tube) Calculations (핀이 있는 상태의 동측 계산) ---
    Af = 2 * Hf * N
    Cs = math.pi * (Dis**2 - Dot**2) / 4
    Nf = Cs - w * Af / 2
    
    # Equivalent diameter for finned section
    Deq_f = 4 * Nf / (math.pi * (Dis + Dot) - N * w + Af)
    Ao = math.pi * Dot + Af
    
    # Fin efficiency calculation (X, e, ep)
    X = Hf * math.sqrt(hifs / (6 * Kt * w))
    # e = tanh(X)/X
    e = (math.exp(X) - math.exp(-X)) / (math.exp(X) + math.exp(-X)) / X
    ep = e * Af / Ao + 1 - Af / Ao
    
    # Finned-side heat transfer coefficient
    hifd = hifs * ep
    
    # --- Overall Heat Transfer Coefficient (총괄 열전달 계수) ---
    Ar = math.pi * Dit / (math.pi * Dot + Af)
    Uo = 1 / ((Dot - Dit) / Kt + 1 / hift * Ar + 1 / hifd)
    
    # --- Output (결과 출력) ---
    print(f"Shell-side heat transfer coefficient: {hifd:10.7f}")
    print(f"Tube-side heat transfer coefficient:  {hift:10.7f}")
    print(f"Fin efficiency:                     {ep:10.7f}")
    print(f"Overall heat transfer coefficient:   {Uo:10.7f}")

if __name__ == "__main__":
    calculate_htcoef()