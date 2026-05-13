import numpy as np

def catLH(y, De, n, h, k, Kr, Ca0):
    """
    Langmuir-Hinshelwood 동역학을 따르는 촉매 펠릿 내 확산 및 반응 모델
    y: 각 격자점에서의 무차원 농도 배열
    De: 유효 확산 계수
    n: 격자 분할 수
    h: 격자 간격
    k: 반응 속도 상수
    Kr: 흡착 상수
    Ca0: 표면 농도 (경계 조건)
    """
    m = n + 1
    # 결과를 저장할 배열 초기화
    z = np.zeros(m)
    # 반경 방향 거리 x 계산
    x = np.array([h * i for i in range(m)])
    
    for i in range(m):
        # Langmuir-Hinshelwood 반응 속도식 계산
        # rxn = k * y / sqrt(1 + Kr * y^2)
        rxn = k * y[i] / np.sqrt(1 + Kr * y[i]**2)
        
        if i == 0:
            # 중심부 경계 조건 (i=1): 대칭성 적용 (dy/dx = 0)
            # z = 2*De*(y[i+1] - y[i])/h^2 - rxn
            z[i] = 2 * De * (y[i+1] - y[i]) / (h**2) - rxn
            
        elif i == m - 1:
            # 표면 경계 조건 (i=m): y = Ca0
            z[i] = y[i] - Ca0
            
        else:
            # 내부 격자점 (i=2 ~ n): 확산 방정식의 차분 형태
            # z = De*(d2y/dx2 + (1/x)*dy/dx) - rxn
            term1 = De * (y[i+1] - 2 * y[i] + y[i-1]) / (h**2)
            term2 = De * (y[i+1] - y[i-1]) / (h * x[i])
            z[i] = term1 + term2 - rxn
            
    return z