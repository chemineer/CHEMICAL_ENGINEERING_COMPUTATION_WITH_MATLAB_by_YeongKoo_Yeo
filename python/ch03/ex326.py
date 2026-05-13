from hform import hform

# 온도 설정 (Kelvin 단위로 가정)
T = 500

# 메탄의 표준 생성열 계산
hf = hform(T, 'methane')

# 결과 출력
print(f"온도 T = {T}K 일 때, 메탄의 표준 생성열(hf) = {hf}")