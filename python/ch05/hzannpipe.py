import numpy as np
import matplotlib.pyplot as plt

def hzannpipe(L, R1, R2, mu, delP):
    """
    수평 환형관(annulus) 내 뉴턴 층류 흐름의 평균 속도 계산
    """
    # 1. 평균 유속 계산 (해석적 해)
    # log 항을 별도로 정의하여 연산 오류 방지
    ln_ratio = np.log(R2 / R1)
    
    # MATLAB 공식과 동일하게 재구성
    # avgv = delP * (R1^2 + R2^2 - (R2^2 - R1^2) / ln(R2/R1)) / (8 * mu * L)
    numerator = delP * (R1**2 + R2**2 - (R2**2 - R1**2) / ln_ratio)
    denominator = 8 * mu * L
    avgv = numerator / denominator
    
    # 2. 분포 계산 및 시각화 (기존 로직 유지)
    n = 100
    r = np.linspace(R1, R2, n)
    h = (R2 - R1) / (n - 1) # h 수정: 간격은 n-1로 나누는 것이 정확함
    
    v = (delP / (4 * mu * L)) * (R2**2 - r**2 + (R2**2 - R1**2) / ln_ratio * np.log(r / R2))
    
    # 전단 응력 분포
    dv = np.gradient(v, h) # np.diff 대신 np.gradient 사용 (중앙 차분법으로 더 정확함)
    taurx = -mu * dv
    
    # 시각화
    plt.figure(figsize=(10, 4))
    
    # 속도 분포 그래프
    plt.subplot(1, 2, 1)
    plt.plot(v, r)
    plt.grid(True)
    plt.xlim(0, 0.2)
    plt.ylim(R1, R2)
    plt.xlabel('v_x(r)')
    plt.ylabel('r')
    plt.title('Velocity Distribution')
    
    # 전단 응력 분포 그래프
    plt.subplot(1, 2, 2)
    plt.plot(taurx, r)
    plt.grid(True)
    plt.xlabel('\\tau_{rx}(r)')
    plt.ylabel('r')
    plt.title('Shear Stress Distribution')
    
    plt.tight_layout()
    plt.show()
    
    return avgv

# 사용 예시:
# avg_vel = hzannpipe(10, 0.05, 0.1, 0.001, 500)
# print(f"Average Velocity: {avg_vel:.4f} m/s")