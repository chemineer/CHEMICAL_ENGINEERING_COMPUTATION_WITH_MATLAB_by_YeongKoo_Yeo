from delH import delH

# 매틀랩 코드의 입력값 할당
# 매틀랩의 배열 [1.702 9.081e-3 -2.164e-6 0]을 리스트로 정의
params = [1.702, 9.081e-3, -2.164e-6, 0]
T1 = 533.15
T2 = 873.15

# 함수 호출 (Q와 mc를 각각 반환받음)
Q, mc = delH(params, T1, T2)

# 결과 출력
print(f"엔탈피 변화 (Q) = {Q}")
print(f"평균 열용량 (mc) = {mc}")