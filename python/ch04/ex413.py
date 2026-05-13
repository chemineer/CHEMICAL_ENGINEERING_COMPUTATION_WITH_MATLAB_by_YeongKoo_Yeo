import numpy as np
from scipy.optimize import fsolve

def wilson_parameters():
    # 1. 입력 데이터 설정
    w1 = 78       # 조성 (중량 %)
    mw1 = 46.07   # 성분 1의 분자량
    mw2 = 114.0   # 성분 2의 분자량
    Ta = 77       # 온도 (°C)
    P = 760       # 압력 (mmHg)

    # Antoine 방정식 상수
    A1, B1, C1 = 8.04494, 1554.3, 222.65
    A2, B2, C2 = 6.92374, 1355.126, 209.517

    # 2. 몰 분율 계산
    x1 = (w1 / mw1) / (w1 / mw1 + (100 - w1) / mw2)
    x2 = 1 - x1

    # 3. 증기압 및 활동도 계수 계산
    P1v = 10**(A1 - B1 / (Ta + C1))
    P2v = 10**(A2 - B2 / (Ta + C2))
    gam1 = P / P1v
    gam2 = P / P2v

    # 4. Wilson 식 비선형 방정식 정의
    def wilact(g, x1, x2, gam1, gam2):
        g12, g21 = g
        t1 = x1 + g12 * x2
        t2 = x2 + g21 * x1
        
        # Wilson 식에 따른 잔차 방정식
        f1 = np.log(gam1) + np.log(t1) - (g12 * t2 - g21 * t1) * x2 / (t1 * t2)
        f2 = np.log(gam2) + np.log(t2) + (g12 * t2 - g21 * t1) * x1 / (t1 * t2)
        return [f1, f2]

    # 5. fsolve를 이용한 매개변수 산출
    g0 = [0.1, 0.1]  # 초기 추정값
    g12, g21 = fsolve(wilact, g0, args=(x1, x2, gam1, gam2))

    print(f"Wilson 매개변수 g12: {g12:.6f}")
    print(f"Wilson 매개변수 g21: {g21:.6f}")

if __name__ == "__main__":
    wilson_parameters()