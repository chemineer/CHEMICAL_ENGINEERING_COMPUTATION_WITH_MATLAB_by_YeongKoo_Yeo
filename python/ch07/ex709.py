import numpy as np
from sinevap import sinevap  # sinevap.py 모듈을 임포트

# 1. 데이터 입력을 위한 클래스 또는 딕셔너리 정의
# MATLAB의 구조체(struct)는 파이썬의 딕셔너리(dict)로 변환하는 것이 가장 일반적입니다.
evdat = {
    'mf': 29000,    # 공급 액체 질량 유량
    'Tf': 60,       # 공급 온 도
    'xf': 0.25,     # 공급 농도
    'xp': 0.6,      # 제품 농도
    'Ps': 25,       # 증기 압력
    'Pv': 1.69,     # 증발기 내부 압력
    'U': 300        # 총괄 열전달 계수
}

# 2. sinevap 함수 호출 및 결과 출력
# sinevap.py 내에 sinevap(evdat) 함수가 정의되어 있어야 합니다.
try:
    res = sinevap(evdat)
    
    # 결과 출력 (MATLAB의 res 결과 확인)
    print("--- Evaporation Design Results ---")
    print(res)
    
except AttributeError:
    print("Error: sinevap.py 내에 'sinevap' 함수가 정의되어 있지 않습니다.")
except Exception as e:
    print(f"An error occurred: {e}")