import numpy as np
import multievapSI  # multievapSI.py 모듈 임포트

# 1. 입력 데이터 설정 (MATLAB의 struct를 dict로 변환)
# 각 변수명과 값은 제공된 소스 코드의 설정을 그대로 따릅니다.
evdat = {
    'xp': 0.60,                                  # 제품 농도
    'xf': 0.15,                                  # 공급 농도
    'Ps': 205602.9,                              # 증기 압력
    'Pn': 8756,                                  # 마지막 효용 증발기 압력
    'mf': 20412,                                 # 공급 액체 질량 유량
    'Tf': 15.6,                                  # 공급 온도
    'U': np.array([10.834, 7.155, 3.986]) * 1e6  # 각 효용별 총괄 열전달 계수
}

# 2. multievapSI 함수 호출
# multievapSI.py 내에 multievapSI(evdat) 함수가 정의되어 있어야 합니다.
try:
    results = multievapSI.multievapSI(evdat)
    print(results)
    
    # 3. 결과 출력 (results.mv 값 확인)
    if hasattr(results, 'mv'):
        # results가 객체(class instance)인 경우
        print(f"Total Vapor Flow (mv): {results.mv}")
    elif isinstance(results, dict) and 'mv' in results:
        # results가 딕셔너리(dict)인 경우
        print(f"Total Vapor Flow (mv): {results['mv']}")
    else:
        print("결과에 'mv' 항목이 존재하지 않습니다.")
        print("전체 결과:", results)

except AttributeError:
    print("Error: 함수가 정의되어 있지 않습니다.")
except Exception as e:
    print(f"오류 발생: {e}")