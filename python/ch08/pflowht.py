import numpy as np

def pflowht(x, T, pf):
    """
    파이프를 흐르는 유체의 열전달을 계산합니다.
    
    Inputs:
      x: 위치 벡터 (함수 서명에 포함되어 있으나 본문에서 사용되지 않음)
      T: 현재 온도 상태 벡터 (n,)
      pf: 매개변수 객체 (r, v, n, h, alpa, Tb 속성 포함)
      
    Outputs:
      dT: 온도 변화율 벡터 (n,)
    """
    # 데이터 추출
    r = pf.r
    v = pf.v
    n = pf.n
    h = pf.h
    alpa = pf.alpa
    Tb = pf.Tb
    
    dT = np.zeros(n)
    
    # 차분 방정식 모델
    for k in range(n):  # MATLAB의 1:n 루프를 0:n-1로 변환
        # MATLAB 인덱스 1 기반을 0 기반으로 매칭
        # k=0 (MATLAB의 k=1)
        if k == 0:
            s = 2 * (T[k+1] - T[k]) / h**2
            d = 0
        # k=n-1 (MATLAB의 k=n)
        elif k == n - 1:
            s = (Tb - 2 * T[k] + T[k-1]) / h**2
            d = (Tb - T[k-1]) / (2 * h * r[k+1]) # r 인덱스 주의: r이 n+1 크기라고 가정
        # 그 외 내부 점
        else:
            s = (T[k+1] - 2 * T[k] + T[k-1]) / h**2
            d = (T[k+1] - T[k-1]) / (2 * h * r[k])
            
        dT[k] = (alpa / v[k]) * (s + d)
        
    return dT