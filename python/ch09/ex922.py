import control as ct
import matplotlib.pyplot as plt

# 1. 전달함수의 분자(num)와 분모(den) 정의
# 매트랩: num = [1 3]; den = [1 2 0];
num = [1, 3]    # s + 3
den = [1, 2, 0] # s^2 + 2s

# 2. 시스템(전달함수) 생성
sys = ct.tf(num, den)

# 3. 근궤적(Root Locus) 그리기
# ct.root_locus는 근궤적 데이터를 계산하고 그래프를 생성합니다.
plt.figure(figsize=(8, 6))
ct.root_locus(sys, grid=True)

# 4. 그래프 설정 (매트랩의 title, xlabel, ylabel 대응)
plt.title('Root locus')
plt.xlabel('Real axis')
plt.ylabel('Imaginary axis')

# 그래프 출력
plt.show()