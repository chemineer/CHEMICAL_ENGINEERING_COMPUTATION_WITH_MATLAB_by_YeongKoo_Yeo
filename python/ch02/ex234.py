from scipy.interpolate import interp1d, PchipInterpolator, CubicSpline

T_cp = np.arange(373, 874, 100)
Cp_data = np.array([29.189, 29.291, 29.462, 29.678, 29.971, 30.269])
target_T = 580

lCp = interp1d(T_cp, Cp_data, kind='linear')(target_T)
pCp = PchipInterpolator(T_cp, Cp_data)(target_T)
sCp = CubicSpline(T_cp, Cp_data)(target_T)
nCp = interp1d(T_cp, Cp_data, kind='nearest')(target_T)

print(f"linear: {lCp:.4f}\npchip: {pCp:.4f}\nspline: {sCp:.4f}\nnearest: {nCp:.4f}")