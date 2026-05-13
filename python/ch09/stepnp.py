import numpy as np
from scipy import signal
from math import factorial

def stepnp(num, den, t0, delt, fint, ms):
    """
    SISO 시스템의 계단 응답을 계산합니다 (MATLAB 원본 로직 이식).
    
    Parameters:
    num  : 전달 함수의 분자 계수
    den  : 전달 함수의 분모 계수
    t0   : 단위 계단 입력이 도입되는 시간
    delt : 시간 간격 (Time step)
    fint : 최종 응답 시간
    ms   : 계단 입력의 크기 (Step size)
    
    Returns:
    y : 계단 응답 결과 배열
    t : 시간 배열
    """
    
    # 1. (전달 함수) * (계단 입력 1/s)의 부분 분수 분해
    # MATLAB의 conv(den, [1 0])는 분모에 s를 곱하는 것과 같습니다.
    den_step = np.convolve(den, [1, 0])
    r, p, k = signal.residue(num, den_step)
    
    # 2. 계산 시간 간격 설정
    t = np.arange(t0, fint + delt, delt)
    
    # 3. 극점의 중복도(Multiplicity) 식별
    # MATLAB 코드의 로직을 파이썬 방식으로 구현
    p_len = len(p)
    mult = np.zeros(p_len, dtype=int)
    for j in range(p_len):
        n = 0
        for i in range(p_len):
            # 부동 소수점 오차를 고려하여 근사값으로 비교
            if np.isclose(p[j], p[i]):
                n += 1
        mult[j] = n

    # 4. 역 라플라스 변환을 이용한 계단 응답 계산
    y = np.zeros(len(t), dtype=complex) # 복소수 극점 대응을 위해 complex 타입 사용
    j = 0
    while j < p_len:
        m = mult[j]
        for i in range(1, m + 1):
            # MATLAB: r(j+i-1) * ((t-t0)^(i-1)) * exp(p(j)*(t-t0)) / (i-1)!
            term = r[j + i - 1] * ((t - t0)**(i - 1)) * np.exp(p[j] * (t - t0)) / factorial(i - 1)
            y += term
        j += m # 중복도만큼 인덱스 건너뛰기
    
    # 입력 크기(ms) 반영 및 실수부만 취함 (허수부는 계산 오차임)
    y = ms * np.real(y)
    
    return y, t

# --- 사용 예시 ---
if __name__ == "__main__":
    # 예: G(s) = 1 / (s^2 + 3s + 2)
    num = [1]
    den = [1, 3, 2]
    t0, delt, fint, ms = 0, 0.1, 10, 1.0
    
    y, t = stepnp(num, den, t0, delt, fint, ms)
    
    print(f"Time (first 5): {t[:5]}")
    print(f"Response (first 5): {y[:5]}")