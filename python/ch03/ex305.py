from denL import denL
import compID # denL 내부에서 compID를 사용할 경우를 대비해 임포트

# 1. 액체 밀도 계산 함수 호출
# MATLAB: rw = denL(150, 'Water')
# 온도 150도에서 물(Water)의 밀도를 계산합니다.
rw = denL(150, 'Water')

# 2. 결과 출력
print(f"Density of Water at 150: {rw}")