import numpy as np

def mArea(x, D, E, F, R, S, T, alpa, r, xf, uf):
    """
    MATLAB mArea 함수를 파이썬으로 변환한 코드입니다.
    입력값 x, D, E 등이 NumPy 배열이거나 스칼라일 때 모두 작동합니다.
    """
    
    # 공통으로 사용되는 제곱근 항 계산
    common_sqrt = np.sqrt((D**2) * (x**2) + 2 * E * x + F**2)
    
    # u 계산
    u = -D * x + common_sqrt
    
    # fi 계산
    fi = (D * x - F) + common_sqrt
    
    # ud 계산 (MATLAB의 .^ 연산은 파이썬의 ** 연산에 대응)
    term1 = (1 - xf)
    term2 = ((uf - E/D) / (u - E/D))**R
    term3 = ((uf - alpa + F) / (u - alpa + F))**S
    term4 = ((uf - F) / (u - F))**T
    
    ud = term1 * term2 * term3 * term4
    
    # 최종 결과 y 계산
    # MATLAB의 ./ 연산은 NumPy에서 배열 간 / 연산으로 처리됩니다.
    denominator = (fi - x) * (1 / (1 + x) - r / (1 + fi))
    y = ud / denominator
    
    return y