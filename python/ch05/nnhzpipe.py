import numpy as np
import matplotlib.pyplot as plt

def nnhzpipe(L, R, K, n, delP):
    """
    비뉴턴 유체의 수평 파이프 흐름에 대한 평균 속도 및 속도 프로파일 계산
    
    input:
    L: 파이프 길이 (m)
    R: 파이프 반경 (m)  # MATLAB 원본 주석에 직경으로 되어있으나, 코드 로직상 반경임
    K: 유체 일관성 지수 (N*s^n/m^2)
    n: 유동 지수
    delP: 압력 강하 (Pa)
    """
    
    # h: 단계 크기, r: 반경 위치 배열
    h = 2 * R / 100
    r = np.arange(-R, R + h, h)
    
    # 속도 프로파일 및 평균 속도 계산
    Rn = R ** ((n + 1) / n)
    Pn = (delP / (2 * K * L)) ** (1 / n)
    
    v = (1 - (np.abs(r) / R) ** ((n + 1) / n)) * Pn * Rn * n / (n + 1)
    avgv = Rn * Pn * n / (3 * n + 1)
    
    # 전단 응력 계산
    dvr = np.diff(v) / h
    dvr = np.append(dvr, dvr[-1])
    taurx = -K * np.abs(dvr) ** (n - 1) * dvr
    
    # 그래프 시각화
    plt.figure(figsize=(10, 4))
    
    plt.subplot(1, 2, 1)
    plt.plot(v, r)
    plt.grid(True)
    plt.xlabel('v_x(r)')
    plt.ylabel('r')
    plt.axis('tight')
    
    plt.subplot(1, 2, 2)
    plt.plot(taurx, r)
    plt.grid(True)
    plt.xlabel('\\tau_{rx}(r)')
    plt.ylabel('r')
    plt.axis('tight')
    
    plt.tight_layout()
    plt.show()
    
    return avgv, taurx