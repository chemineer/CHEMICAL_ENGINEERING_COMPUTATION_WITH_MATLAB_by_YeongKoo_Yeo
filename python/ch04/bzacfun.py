import numpy as np

def bzacfun(x, y):
    # solve_ivp에서 y는 1차원 배열(array)로 들어오므로 스칼라 값으로 추출
    y_val = y[0] if isinstance(y, np.ndarray) else y
    
    # MATLAB 수식과 완벽히 동일하게 적용 (불필요한 if문 제거)
    numerator = (y_val * (1 - y_val) / (y_val - x)) * (
        -1463.0572 * x**3 + 2745.4788 * x**2 - 1788.1398 * x + 556.8470
    )
    
    denominator = (
        -365.7643 * x**4 + 915.1596 * x**3 - 894.0699 * x**2 + 556.8470 * x + 56.2570
    )
    
    dy = numerator / denominator
    
    # solve_ivp 규격에 맞게 리스트로 반환
    return [dy]