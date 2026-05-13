import numpy as np
from scipy.optimize import fmin
from types import SimpleNamespace
from rfobj import rfobj  # 작성하신 rfobj.py 임포트

# 1. 데이터 설정 (opdat 객체 생성)
opdat = SimpleNamespace(
    we=0.8, eta=0.75, S=12000, K=0.05,
    rhoG=2, rhoL=850, rhos=8000,
    lamb=800, lambs=1800, Css=0.05, Cst=10,
    F=100, T=70, P=760, alpa=2.3,
    xB=0.05, xD=0.85, xF=0.4
)

# 2. 파라미터 추출
alpa = opdat.alpa
xD = opdat.xD
xF = opdat.xF

# 3. 최적 환류비(rfopt) 찾기
rf0 = 1.5  # 초기 추정값
# fmin(fminsearch 대응)을 사용하여 rfobj 함수의 최솟값을 찾음
# rfobj의 첫 번째 인자 외의 추가 인자(opdat)는 args 파라미터로 전달
rfopt_res = fmin(rfobj, rf0, args=(opdat,), disp=False)
rfopt = rfopt_res[0]

# 4. Fenske 식에 의한 최소 환류비 및 경험칙(Rule of thumb) 계산
minrf = (xD/xF - alpa*(1-xD)/(1-xF)) / (alpa-1)
rfpr = 1.2 * minrf

# 5. 결과 출력
print(f"Optimal reflux ratio = {rfopt:.6g}")
print(f"Reflux ratio by the rule of thumb = {rfpr:.6g}")