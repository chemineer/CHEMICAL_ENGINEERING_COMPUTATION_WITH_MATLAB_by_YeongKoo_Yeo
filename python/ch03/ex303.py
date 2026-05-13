from satH2Oprop import satH2Oprop

# 1. 온도(T=150)를 기준으로 포화 증기('V')의 물성치 계산
# MATLAB: x = satH2Oprop('T','V',150)
x_vapor = satH2Oprop('T', 'V', 150)
print("Saturated Vapor Properties at T=150:")
print(x_vapor)

# 2. 온도(T=150)를 기준으로 포화 액체('L')의 물성치 계산
# MATLAB: x = satH2Oprop('T','L',150)
x_liquid = satH2Oprop('T', 'L', 150)
print("\nSaturated Liquid Properties at T=150:")
print(x_liquid)