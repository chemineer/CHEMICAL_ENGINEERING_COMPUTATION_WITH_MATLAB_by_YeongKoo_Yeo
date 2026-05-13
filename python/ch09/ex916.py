import matplotlib.pyplot as plt
from secpro import secpro # secpro.py에서 함수 임포트

# 1. 그래프 생성 준비 (매트랩의 hold on 효과를 위함)
plt.figure(figsize=(8, 6))

# 2. Kc 값을 변경하며 secpro 함수 호출
# 매트랩 코드의 의도에 따라 Kc 값을 0.5, 1, 2로 설정하여 실행합니다.
kc_values = [0.5, 1, 2]

for kc in kc_values:
    # secpro.py 내의 secpro 함수를 실행하여 그래프를 겹쳐 그림
    secpro(kc)

# 3. 특정 위치에 텍스트 표시 (매트랩의 text 함수 대응)
plt.text(2.1, 1.1, 'Kc=0.5')
plt.text(4.1, 0.86, 'Kc=1')
plt.text(5.1, 0.68, 'Kc=2')

# 4. 최종 그래프 출력
plt.title('Step Responses of 2nd-order Process (ex916)')
plt.show()