import numpy as np

def pitzersg(w, Pc, Tc, T):
    """
    Estimation of surface tension using Pitzer's relation
    w: acentric factor
    Pc: critical pressure (bar)
    Tc: critical temperature (K)
    T: temperature (K)
    output: surface tension (dyne/cm)
    """
    
    # 환산 온도 계산
    Tr = T / Tc
    
    # Pitzer's relation 수식 적용
    # MATLAB의 요소별 연산에 대응하기 위해 numpy 배열 연산 사용
    term1 = Pc ** (2/3)
    term2 = Tc ** (1/3)
    term3 = (1.86 + 1.18 * w) / 19.05
    term4 = ((3.75 + 0.91 * w) / (0.291 - 0.08 * w)) ** (2/3)
    term5 = (1 - Tr) ** (11/9)
    
    sg = term1 * term2 * term3 * term4 * term5
    
    # 결과 출력
    print(f'Surface tension = {sg:g} dyne/cm')
    
    return sg