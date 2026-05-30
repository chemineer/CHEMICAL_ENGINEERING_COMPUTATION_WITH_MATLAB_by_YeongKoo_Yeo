## 9 Process Control

### Example 9.1 Laplace Transform 
Find the Laplace transform of 
$f(t) = 1 + t+ t^2+sin at -t cos bt $. 
```py
import sympy as sp

# 심볼 정의
t, a, b = sp.symbols('t a b')

# 함수 정의
f = 1 + t + t**2 + sp.sin(a*t) - t*sp.cos(b*t)

# 라플라스 변환 수행
Lf = sp.laplace_transform(f, t, sp.symbols('s'))

# 결과 출력 (라플라스 변환된 식과 수렴 조건이 함께 나옵니다)
print(Lf[0])
```

### Example 9.2 Inverse Laplace Transform
Find the inverse Laplace transform of 
$F(s ) = 2s /(s^2 + 4s +1) $

### Example 9.3 Partial Fraction Expansion
$F(s )=(s^3 +5s^2 +9s +7)/(s^2 +3s +2)$ can be expanded into 
$F(s )= -1/(s +2)+2/(s +1)+(s +2) $. Verify this using MATLAB.

### Example 9.4 Representation of the Transfer Function
Use MATLAB to display $G(s )=(2s +1)/(s^2 +3s +2)=(2s +1)/((s +1)(s +2))$. 

### Example 9.5 Use of the Built-In Function tf
Find the transfer function G when the numerator is given by 2s +1 and the denominator is
given by $s^2 +3s +2 $

### Example 9.6 Overall Transfer Function
Find the overall transfer function for the simple feedback loop shown in Figure 9.2. In 
Figure 9.2, $G =(2s +1)/(s^2 +3s +2)$, H=1/(s +1).
![그림](f902.png)

### Example 9.7 From Transfer Function to State-Space Model
Find a state-space model for the process whose transfer function is given by
$$
G=\frac{2s+1}{s^3+3s^2+2s+1}
$$

### Example 9.8 State-Space Representation 
Find a state-space representation for the system given by
$$
\frac{d^2y}{dt^2}+1.5\frac{dy}{dt}+y=\frac{du}{dt}+2u, y(0)=\frac{dy}{dt}(0)=u(0)=0
$$

### Example 9.9 Step Response of a 1st-Order Process
Plot the unit step response of a 1st-order process during the time interval [0, 10]. The transfer function of the process is given by 3/(2s + 1)

### Example 9.10 Sinusoidal Response of a 1st-Order Process
Plot the response curve of the 1st-order system 3/(2s + 1) to the sinusoidal input u = sin(3t )
during the time interval [0, 10]. 

### Example 9.11 Dynamics of Two Heated Tanks in Series
Figure 9.7 shows two heated tanks in series. The two tanks are connected with pipes at the bottom and at a certain height H. If the liquid level of the first tank h1 is greater than H, the liquid flows through both connecting pipes. If h1 < H , the liquid flows through the bottom pipe only. Heats Q1 and Q2 are introduced into tanks from external heat sources. The cross-sectional areas of tank 1 and tank 2 are A1 and A2, respectively. From the material and energy balances for this system, we obtain 

$$
A_1 \frac{dh_1}{dt} = F_0 -F_{1t} -F_{1b}, A_2 \frac{dh_2}{dt} = F_{1t} +F_{1b} -F_{2}
$$
$$
\frac{dT_1}{dt} =\frac{F_0}{A_1h_1}(T_0-T_1) +\frac{Q_1}{\rho C_pA_1h_1}, \frac{dT_2}{dt} =(F_{1t} +F_{1b})(T_1-T_2) +\frac{Q_2}{\rho C_pA_2h_2}
$$

The flow rates in these equations are given by 
$$
F_{1b}=c_1\sqrt{h_1-h_2}, F_{2}=c_2\sqrt{h_2}, F_{1t}=\left\{\begin{matrix}
0 : h_1 < H \\
c_1 \sqrt{h_1-H} : h_1 > H, h_2 < H \\
c_2 \sqrt{h_1-h_2} : h_2 > H \end{matrix}\right.
$$

where c1 and c2 are valve coefficients. Initially, both liquid tanks contain  100 l of water at 25°C. 
Generate time profiles of temperatures (T1, T2) and liquid levels (h1, h2) for 0 < t < 10(min )using the given data. 

Data: 
A1 = A2 = 0.25 m2,
F0 = 0.4 m3/min
T0 = 25 °C
c1 = c2 = 0.6 m2.5/min
ρCp = 4180 kJ/kg,
Q1 = Q2 = 6000 kJ/min, H = 0.5 m. 

### Example 9.12 Step Response of a 2nd-Order Process 
In the 2nd-order process represented by $1/(\tau ^2s^2 +2\tau \xi s +1)  $
, the value of the time constant is τ = 0.5
. Plot the step response curves when ξ is 0.5, 1.0, and 1.5 during the time interval [0, 10].

### Example 9.13 Step Response of a Higher-Order Process
Plot the unit step response curves for the higher-order process represented by $1/(2s +1)^2 $, $1/(2s +1)^4 $, $1/(2s +1)^5 $ during the time interval [0, 20].

### Example 9.14 Step Response of a 1st-Order Plus Time Delay Process
Plot the unit step response curve of a 1st-order process with time delay, the transfer function of which is given by 
$$
G_1(s)=\frac{Y(s)}{X(s)}=\frac{3exp(-1.6s)}{3s+1}
$$
Compare the result with the step response of the 4th-order process given by
$$
G_2(s)=\frac{Y(s)}{X(s)}=\frac{3}{(0.1s+1)(0.5s+1)(s+1)(3s+1)}
$$

### Example 9.15 Step Response of a Feedback Control System
The overall closed-loop transfer function of a feedback control system is given by 
$$
C=\frac{K_c}{5s+1+K_c}R
$$
Calculate and plot the closed-loop response to a unit step change in the set point for three values of the proportional controller gain: Kc = 5, 20, and 50.

### Example 9.16 Step Responses for Proportional Control of a 2nd-Order Process
A proportional controller is to be used to control a 2nd-order process given by G= 0.5/( s(0.5s + 1)). Calculate and plot the closed-loop response to a unit step change in the set point for three values of the proportional controller gain: Kc = 0.5, 1 and 2. Assume that Kv = Km= =1. 

### Example 9.17 Step Response of a Feedback Control System Using a Proportional-Integral Controller 
A proportional-integral (PI) controller is used in a feedback control system. The process transfer function is given by 
Gp(s) = 5/((s+1)(2s+1))
, the gain of the control valve is Kv = 0.01, and the gain of the sensor/transducer is Km = 20. Use Simulink to find the step response curve to a unit step change in the set point. The transfer function of the PI controller is given by Gc (s) = 2(1 + 1/(5s )). 


### Example 9.18 Proportional Integral Control of a Batch Reactor
The following exothermic consecutive reactions are carried out in a batch reactor fitted with a cooling coil through which cooling water is passed to remove the exothermic heat, as shown in Figure 9.19: 

$\ce{A ->[k1] B ->[k2] C} $
From the material balances for species A and B, we obtain 
$dC_A/dt=-k_1C_A^2, dC_B/dt=k_1C_A^2-k_2C_B $
where 
k1 and k2 are the reaction rate constants
$C_A$ and $C_B$ are the concentrations of species A and B, respectively 
k1 and k2 are represented as 
$k1=A_1exp(-E_1/(RT)), k2=A_2exp(-E_2/(RT)) $
The energy balance for the batch reactor gives 
$\frac{dT}{dt} =\frac{-\Delta H_1}{\rho C_p}k_1C_A^2 +\frac{-\Delta H_2}{\rho C_p}k_2C_B +\frac{U_j A_j}{\rho C_p V}(T_s-T) -\frac{U_c A_c}{\rho C_p V}(T-T_c) $
where 
($-\Delta H_1$ ) is the heat of reaction for A -> B
($-\Delta H_2$ ) is the heat of reaction for B -> C
Ts and Tc are the steam and coolant temperatures, respectively 
Uj and Uc are the overall heat transfer coefficients of the jacket and coolant, respectively

Uc is assumed to be a function of the coolant flow rate Fc as 
$\frac{1}{U_c} =\frac{1}{4550 F_c^{0.8}} + \frac{1}{10.8} $

For the present case, the reactor temperature should precisely followed the desired trajectory 
given by 
$T_d(t)=54+71 exp(-0.0025 t) $

One way to control T(t ) is to introduce a parameter u defined by 
$T_s=(T_{s,max}-T_{s,min})u +T_{s,min}, U_c=(U_{c,max}-U_{c,min})u +U_{c,max} $

u=0 represents the maximum cooling and u=1 the maximum heating of the system. Substituting 
these relations into the energy balance, followed by rearrangement, yields 
$\frac{dT}{dt} =\gamma_1k_1C_A^2 +\gamma_2k_2C_B +(a_1 +A_2T) +(b_1 +b_2 T)u $

where 
$\gamma_1 =\frac{(-\Delta H_1)}{\rho C_p}, \gamma_2 =\frac{(-\Delta H_2)}{\rho C_p}, a_1 =\frac{U_jA_jT_{s,min}+U_{c,max}A_cT_c}{\rho C_p V}, a_2 =-\frac{U_jA_j+U_{c,max}A_c}{\rho C_p V} $

$b_1 =\frac{U_jA_j(T_{s,max}-T_{s,min})-(U_{c,max}-U_{c,min})A_cT_c}{\rho C_p V}, b_2 =\frac{(U_{c,max}-U_{c,min})A_c}{\rho C_p V} $

A simple PI controller is used in the control, and uis determined by 
$u(t) =u_s +K_c(e(t)+\frac{1}{\tau_I}\int_{0}^{t}e(t)dt ) $

where the control error e(t ) is given by e(t )= Td(t ) -T(t ) and the controller gain is arbitrarily chosen as $K_c =0.1(°C^{-1})$ and τI =360(sec ). Plot CA, CB, Fc, Td, Ts, and T as a function of reaction time t (0 < t < 4000 ). 
Data: CA0  =1.0 kmol/m3 , CB0  =0.0 kmol/m3, A1 =1.1 m3/(kmol sec ) , A2  =172.2 sec -1, E1=
 2.09×10 4 kJ /kmol , E2 =4.18×10 4 kJ /kmol , (-ΔH1 )= 4.18×10 4 kJ /kmol , (-ΔH2 )=8.36× 10 4  kJ /kmol , ρ =1000 kg/ m3, Tc =25°C , Uj =1.16 KJ/(m2 °C sec) , Ucmax = 4.42 KJ/(m2 °C sec ) , 
Ucmin =1.39 KJ/(m2 °C sec) , Tsmax =150°C , Tsmin =70°C , R =8.314 KJ/(kmol K ), Ac /V = 
17 m2/m3, Aj /V =30 m2/m3, Cp =1.0 KJ/(kg °C) , us =1.0 , $K_c =0.1(°C^{-1})$ and τI =360(sec ) . 

### Example 9.19 Step Response of a Feedback Control System Using a PID Controller 
Figure 9.21 shows a simple feedback control system where the process transfer function is 
$G_p(s) = 3/(4s^2+s+1)   $
and a proportional-integral-derivative (PID) controller is used. The controller transfer function is 
$G_c(s) = K_c(1+1/(\tau_I s)+\tau_D s)   $
. 
Calculate and plot the closed-loop response to a unit step change in the set point for various values of τI and τD while the value of Kc is kept constant. 

### Example 9.20 Control of a Stirred Tank Heating Process 
A continuous-stirred tank heating process consists of a stirred tank, heater, and PI controller. 
The liquid feed with density ρ =980 kg/ m3 and heat capacity of Cp = 1.6 kJ/(kg ∙ °C) flows into the heated tank at a constant mass flow rate of w = 250 kg/min and temperature T=50°C i . The volume of the tank is V =2.5 m3. This stream is to be heated to a higher set-point temperature TR =75°C . 
The outlet temperature is measured by a thermocouple as Tt, and the heat flux Q  (kJ /min ) supplied by the heater is adjusted by a PI controller with K =60 kJ/(min °C) and  τI=1.8 min . The thermocouple exhibits 1st-order dynamics with time constant τt =3 min and unity gain. The time delay between the heating tank and the thermocouple is θ =1.2 min. 
(1) The system is initially operating at steady state at a temperature of 75°C(set point). The inlet temperature Ti is suddenly changed to 30°C at time t  =10 min. Assuming no control actions (open loop, Kc =0 ), plot the profiles of the fluid temperature in the tank T, the measured temperature Tt, and the fluid outlet temperature T1. 
(2) If the controller is engaged (closed loop), plot the profiles of T, Tt, and T1. 
(3) If the PI controller is replaced by a proportional (P) controller with Kc =60 kJ/(min °C) , plot the profiles of T, Tt, and T1.

### Example 9.21 Root Locus
Plot the trajectory of poles when the controller gain Kc changes from 1 to 40 for the system with the characteristic equation given by 
$$
s^3+6s^2+11s+6+2Kc =0
$$

### Example 9.22 Root Locus
Plot the root locus diagram for a system with the characteristic equation given by
$$
1+ \frac{K(s+3)}{s^2+2s} =0
$$

### Example 9.23 Bode Diagram of a 2nd-Order Process 
Plot the Bode diagram of a 2nd-order process 
$$
G(s)=1/(2.25 s^2+3 ξs +1) 
$$
when ξ is 0, 0.25, 0.5, 0.75 and 1.

### Example 9.24 Bode Diagram of a 3rd-Order Process 
Plot the Bode diagram of a 3rd-order process with time delay: 
$$
G(s)=\frac{(0.4s+1)exp(-0.2s)}{(0.3s+1)(S+1)^2}
$$

### Example 9.25 Nyquist Diagram of a 2nd-Order Process 
Plot the Nyquist diagram for a 2nd-order process whose transfer function is 
$$
G(s)=\frac{2}{(10s+1)(2.5s+1)}
$$

### Example 9.26 Nyquist Diagram of a Time Delay 
Plot the Nyquist diagram for a 1st-order process with a time delay whose transfer function is
$$
G(s)=\frac{12.76 exp(-s)}{5s+1}
$$

### Example 9.27 Nichols Chart
Plot the Nichols chart for a 2nd-order process whose transfer function is
$$
G(s)=\frac{2}{(10s+1)(2.5s+1)}
$$


### Example 9.28 Ultimate Gain 
The characteristic equation of a feedback control system using a proportional controller is given by
$$
1+K_c\frac{0.8 exp(-2s)}{5s+1}
$$
where Kc is the controller gain. Find the ultimate gain. 
