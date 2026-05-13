import numpy as np

def pbrmult(W, x, ka, kc, Ct0, Ft0, alpha):
    """
    MATLAB의 pbrmult 함수를 파이썬으로 구현
    W: 가중치(여기서는 독립 변수)
    x: 상태 벡터 (x[0]=Fa, x[1]=Fb, x[2]=Fc, x[3]=Fd, x[4]=y)
    ka, kc, Ct0, Ft0, alpha: 모델 매개변수
    """
    # 총 몰유량 계산
    Ft = x[0] + x[1] + x[2] + x[3]
    
    # 농도(C) 계산 (MATLAB의 for 루프를 리스트 컴프리헨션으로 변환)
    # C[0]=Ca, C[1]=Cb, C[2]=Cc, C[3]=Cd
    C = [Ct0 * x[i] * x[4] / Ft for i in range(4)]
    
    # 미분 방정식 값 계산
    # x[0]' = dFa/dW
    # x[1]' = dFb/dW
    # x[2]' = dFc/dW
    # x[3]' = dFd/dW
    # x[4]' = dy/dW
    
    dFa = -ka * C[0] * C[1]**2 - 2 * kc * C[0]**2 * C[2]**3 / 3
    dFb = -2 * ka * C[0] * C[1]**2
    dFc = ka * C[0] * C[1]**2 - kc * C[0]**2 * C[2]**3
    dFd = kc * C[0]**2 * C[2]**3 / 3
    dy = -alpha * Ft / (2 * x[4] * Ft0)
    
    frx = np.array([dFa, dFb, dFc, dFd, dy])
    
    return frx