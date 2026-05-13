import numpy as np
import matplotlib.pyplot as plt

def hzpipe(L, R, mu, delP):
    """
    수평 원형 파이프 내 층류 흐름의 평균 속도 및 속도 분포를 계산합니다.
    
    입력:
    L: 파이프 길이 (m)
    R: 파이프 반지름 (m)
    mu: 점도 (kg/m/s)
    delP: 압력 강하 (Pa)
    
    출력:
    v: 평균 속도 (m/s)
    """
    # 속도 분포를 위한 r 좌표 생성
    r = np.linspace(-R, R, 100)
    
    # 속도 분포 (vr) 및 평균 속도 (v) 계산
    # vr = (delP * R^2 * (1 - (r/R)^2)) / (4 * mu * L)
    vr = (delP * R**2 * (1 - (r / R)**2)) / (4 * mu * L)
    v = (delP * R**2) / (8 * mu * L)
    
    # 시각화
    plt.figure()
    plt.plot(vr, r)
    plt.grid(True)
    plt.xlabel('v_x(r)')
    plt.ylabel('r')
    plt.xlim(0, 1)  # 원본 코드의 axis([0 1 -R R]) 반영
    plt.ylim(-R, R)
    plt.title('Velocity Distribution in Horizontal Pipe')
    plt.show()
    
    return v

# 사용 예시:
# avg_v = hzpipe(10, 0.05, 0.001, 500)
# print(f"Average Velocity: {avg_v:.4f} m/s")