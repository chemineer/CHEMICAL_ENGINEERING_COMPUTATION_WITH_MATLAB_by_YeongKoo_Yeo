from satH2Oprop import satH2Oprop

# 1. 압력(P)을 기준으로 포화 액체('L')의 물성치 계산
# MATLAB: x = satH2Oprop('P', 'L', 10540)
x1 = satH2Oprop('P', 'L', 10540)
print(f"Saturated Liquid Properties at P=10540: {x1}")

# 2. 온도(T)를 기준으로 포화 증기('V')의 물성치 계산
# MATLAB: x = satH2Oprop('T', 'V', 198)
x2 = satH2Oprop('T', 'V', 198)
print(f"Saturated Vapor Properties at T=198: {x2}")