import math

def humidest(Hr, Td):
    """
    Hr: 상대 습도 (Relative humidity, %)
    Td: 건구 온도 (Dry bulb temperature, °C)
    """
    
    # 임계 상수 (Critical properties)
    Tc = 647.096  # 임계 온도 (K)
    Pc = 22064000  # 임계 압력 (Pa)
    
    # 상수 (Constants)
    Rw = 461.512244565
    a = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]
    
    # 온도 파라미터 계산
    Td_k = Td + 273.15  # 섭씨를 켈빈으로 변환
    theta = Td_k / Tc
    tau = 1 - theta
    
    # 포화 수증기압(Saturated pressure) 계산
    # 파이썬 리스트 인덱스는 0부터 시작하므로 a[0]~a[5] 사용
    tw = (Tc / Td_k) * (
        a[0] * tau + 
        a[1] * (tau ** 1.5) + 
        a[2] * (tau ** 3) + 
        a[3] * (tau ** 3.5) + 
        a[4] * (tau ** 4) + 
        a[5] * (tau ** 7.5)
    )
    
    Ps = Pc * math.exp(tw)
    
    # 결과 계산
    actual_pres = (Hr / 100) * Ps  # 실제 수증기압 (Actual vapor pressure)
    absolute_ha = actual_pres * 1000 / (Td_k * Rw)  # 절대 습도 (Absolute humidity)
    
    # 결과 출력
    print(f"Actual vapor pressure: {actual_pres:g} Pa")
    print(f"Absolute humidity = {absolute_ha:g} g/m^3")
    
    # 결과를 객체(딕셔너리) 형태로 반환
    return {
        "pres": actual_pres,
        "Ha": absolute_ha
    }

# 사용 예시:
# result = humidest(50, 25)