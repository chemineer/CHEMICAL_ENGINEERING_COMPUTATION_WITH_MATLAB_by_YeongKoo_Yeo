import numpy as np

def pbrmf(w, x, k, fa0, ca0, ka, kb, kc):
    """
    MATLAB의 pbrmf 함수를 파이썬으로 구현
    w: 가중치(여기서는 변수)
    x: 상태 벡터 (입력 배열)
    k, fa0, ca0, ka, kb, kc: 모델 매개변수
    """
    kfc = k * ca0**2 / fa0
    
    # 각 항 계산
    term1 = kfc * (1 - x[0])**2 / (1 + ka * ca0 * (1 - x[0]))
    term2 = kfc * (1 - x[1])**2 / (1 + ka * ca0 * (1 - x[1]) + kc * ca0 * x[1])
    term3 = kfc * (1 - x[2])**2 / (1 + ka * ca0 * (1 - x[2]) + kb * ca0 * (1 - x[2]))**2
    term4 = kfc * (1 - x[3])**2 / (1 + ka * ca0 * (1 - x[3]) + kb * ca0 * (1 - x[3]) + kc * ca0 * x[3])**2
    
    dxdw = np.array([term1, term2, term3, term4])
    
    return dxdw