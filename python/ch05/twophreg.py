import numpy as np

def twophreg(Bx): 
    # 2-phase flow regimes
    # Bx, By: Baker parameters[cite: 3]
    
    # 입력값이 단일 숫자인 경우를 대비해 numpy 배열로 변환 (np.asfarray 미사용)
    Bx = np.array(Bx, dtype=float)
    logBx = np.log(Bx)
    
    C1 = np.exp(9.774459 - 0.6548 * logBx) #[cite: 3]
    C2 = np.exp(8.67964 - 0.1901 * logBx) #[cite: 3]
    C3 = np.exp(11.3976 - 0.6084 * logBx + 0.0779 * logBx**2) #[cite: 3]
    C4 = np.exp(10.7448 - 1.6265 * logBx + 0.2839 * logBx**2) #[cite: 3]
    C5 = np.exp(14.569802 - 1.0173 * logBx) #[cite: 3]
    C6 = np.exp(7.8206 - 0.2189 * logBx) #[cite: 3]
    
    # MATLAB의 [C1' C2' ...]와 유사하게 행 방향으로 합침[cite: 3]
    By = np.array([C1, C2, C3, C4, C5, C6]).T
    return By