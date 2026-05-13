import numpy as np

# 1. 데이터 입력
t = np.array([0, 2.25, 4.50, 6.33, 8.00, 10.25, 12.00, 13.50, 15.60, 17.85, 
              19.60, 27.00, 30.00, 38.00, 41.00, 45.00, 47.00, 57.00, 63.00])

Cb = np.array([0.3335, 0.2965, 0.2660, 0.2450, 0.2255, 0.2050, 0.1910, 0.1794, 0.1632, 
               0.1500, 0.1429, 0.1160, 0.1053, 0.0830, 0.0767, 0.0705, 0.0678, 0.0553, 0.0482])

# 2. 수치 미분 (Numerical Differentiation)
# np.diff는 차분(difference)을 계산합니다.
dCb_raw = np.diff(Cb) / np.diff(t)
# 마지막 원소를 복사하여 길이를 맞춤
dCb = np.append(dCb_raw, dCb_raw[-1])

# 3. 선형 회귀 분석 준비 (A * x = B 형태)
# 반응 속도식: -dCb/dt = k * Cb^n
# 로그 변환: log(-dCb/dt) = n * log(Cb) + log(k)
n_data = len(t)
A = np.column_stack([np.log(Cb), np.ones(n_data)])
B = np.log(-dCb)

# 4. 최소자승법 (Least Squares Method) 계산
# x = (A^T * A)^-1 * A^T * B
x = np.linalg.inv(A.T @ A) @ A.T @ B

# 5. 결과 추출
n_order = x[0]        # 반응 차수 (n)
k_constant = np.exp(x[1])  # 반응 속도 상수 (k)

# 결과 출력
print(f"반응 차수 (n): {n_order:.4f}")
print(f"반응 속도 상수 (k): {k_constant:.4f}")