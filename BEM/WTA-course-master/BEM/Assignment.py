import matlab.engine
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass

        
eng = matlab.engine.start_matlab()


# example input
V = 12.0
omega = 1.267
pitch = 0.0

# call MATLAB
Rx, FN, FT, P = eng.BEM(V, omega, pitch, nargout=4)

# convert to numpy
Rx = np.array(Rx).flatten()
FN = np.array(FN).flatten()
FT = np.array(FT).flatten()



# ============================================================
# STRUCTURAL DATA
# ONLY VALUES NEEDED FOR 2DOF MODEL
# ============================================================

# radius [m]
r = np.array([
    1.50,1.70,2.70,3.70,4.70,5.70,6.70,7.70,8.70,9.70,
    10.70,11.70,12.70,13.70,14.70,15.70,16.70,17.70,
    19.70,21.70,23.70,25.70,27.70,29.70,31.70,33.70,
    35.70,37.70,39.70,41.70,43.70,45.70,47.70,49.70,
    51.70,53.70,55.70,56.70,57.70,58.70,59.20,59.70,
    60.20,60.70,61.20,61.70,62.20,62.70,63.00
])

# mass density [kg/m]
m = np.array([
    678.935,678.935,773.363,740.550,740.042,592.496,
    450.275,424.054,400.638,382.062,399.655,426.321,
    416.820,406.186,381.420,352.822,349.477,346.538,
    339.333,330.004,321.990,313.820,294.734,287.120,
    263.343,253.207,241.666,220.638,200.293,179.404,
    165.094,154.411,138.935,129.555,107.264,98.776,
    90.248,83.001,72.906,68.772,66.264,59.340,
    55.914,52.484,49.114,45.818,41.669,11.453,10.319
])

# flap stiffness [Nm^2]
EI_flap = np.array([
    18110.00E6,18110.00E6,19424.90E6,17455.90E6,
    15287.40E6,10782.40E6,7229.72E6,6309.54E6,
    5528.36E6,4980.06E6,4936.84E6,4691.66E6,
    3949.46E6,3386.52E6,2933.74E6,2568.96E6,
    2388.65E6,2271.99E6,2050.05E6,1828.25E6,
    1588.71E6,1361.93E6,1102.38E6,875.80E6,
    681.30E6,534.72E6,408.90E6,314.54E6,
    238.63E6,175.88E6,126.01E6,107.26E6,
    90.88E6,76.31E6,61.05E6,49.48E6,
    39.36E6,34.67E6,30.41E6,26.52E6,
    23.84E6,19.63E6,16.00E6,12.83E6,
    10.08E6,7.55E6,4.60E6,0.25E6,0.17E6
])

# edge stiffness [Nm^2]
EI_edge = np.array([
    18113.60E6,18113.60E6,19558.60E6,19497.80E6,
    19788.80E6,14858.50E6,10220.60E6,9144.70E6,
    8063.16E6,6884.44E6,7009.18E6,7167.68E6,
    7271.66E6,7081.70E6,6244.53E6,5048.96E6,
    4948.49E6,4808.02E6,4501.40E6,4244.07E6,
    3995.28E6,3750.76E6,3447.14E6,3139.07E6,
    2734.24E6,2554.87E6,2334.03E6,1828.73E6,
    1584.10E6,1323.36E6,1183.68E6,1020.16E6,
    797.81E6,709.61E6,518.19E6,454.87E6,
    395.12E6,353.72E6,304.73E6,281.42E6,
    261.71E6,158.81E6,137.88E6,118.79E6,
    101.63E6,85.07E6,64.26E6,6.61E6,5.01E6
])

# ============================================================
# STRUCTURAL MODEL
# ============================================================

@dataclass
class BladeSection:
    r: float
    m: float
    EI_flap: float
    EI_edge: float

class StructuralModel:

    def __init__(self, r, m, EI_flap, EI_edge):

        self.r = r
        self.m = m
        self.EI_flap = EI_flap
        self.EI_edge = EI_edge

        self.R = np.max(r)

        # damping ratio
        self.zeta = 0.00477465

        self.M = None
        self.K = None
        self.C = None

    # ========================================================
    # MODE SHAPES
    # ========================================================

    def phi_flap(self, r):

        x = r / self.R

        return (
            0.0622*x**2
            + 1.7254*x**3
            - 3.2452*x**4
            + 4.7131*x**5
            - 2.2555*x**6
        )

    def phi_edge(self, r):

        x = r / self.R

        return (
            0.3627*x**2
            + 2.5337*x**3
            - 3.5772*x**4
            + 2.376*x**5
            - 0.6952*x**6
        )

    # ========================================================
    # SECOND DERIVATIVES
    # ========================================================

    def d2phi_flap(self, r):

        x = r / self.R

        d2_dx2 = (
            2*0.0622
            + 6*1.7254*x
            - 12*3.2452*x**2
            + 20*4.7131*x**3
            - 30*2.2555*x**4
        )

        return d2_dx2 / self.R**2

    def d2phi_edge(self, r):

        x = r / self.R

        d2_dx2 = (
            2*0.3627
            + 6*2.5337*x
            - 12*3.5772*x**2
            + 20*2.376*x**3
            - 30*0.6952*x**4
        )

        return d2_dx2 / self.R**2

    # ========================================================
    # MASS MATRIX
    # ========================================================

    def compute_M(self):

        Mf = 0.0
        Me = 0.0

        for i in range(len(self.r)-1):

            dr = self.r[i+1] - self.r[i]

            phif = self.phi_flap(self.r[i])
            phie = self.phi_edge(self.r[i])

            Mf += self.m[i] * phif**2 * dr
            Me += self.m[i] * phie**2 * dr

        self.M = np.array([
            [Mf, 0.0],
            [0.0, Me]
        ])

        return self.M

    # ========================================================
    # STIFFNESS MATRIX
    # ========================================================

    def compute_K(self):

        Kf = 0.0
        Ke = 0.0

        for i in range(len(self.r)-1):

            dr = self.r[i+1] - self.r[i]

            d2phif = self.d2phi_flap(self.r[i])
            d2phie = self.d2phi_edge(self.r[i])

            Kf += self.EI_flap[i] * d2phif**2 * dr
            Ke += self.EI_edge[i] * d2phie**2 * dr

        self.K = np.array([
            [Kf, 0.0],
            [0.0, Ke]
        ])

        return self.K

    # ========================================================
    # DAMPING MATRIX
    # ========================================================

    def compute_C(self):

        if self.M is None:
            self.compute_M()

        if self.K is None:
            self.compute_K()

        Cff = 2 * self.zeta * np.sqrt(
            self.K[0,0] * self.M[0,0]
        )

        Cee = 2 * self.zeta * np.sqrt(
            self.K[1,1] * self.M[1,1]
        )

        self.C = np.array([
            [Cff, 0.0],
            [0.0, Cee]
        ])

        return self.C

    # ========================================================
    # GENERALIZED FORCES
    # ========================================================

    def generalized_forces(self, Rx, FN, FT):

        # interpolate BEM loads onto structural grid
        FN_interp = np.interp(self.r, Rx, FN)
        FT_interp = np.interp(self.r, Rx, FT)

        Qf = 0.0
        Qe = 0.0

        for i in range(len(self.r)-1):

            dr = self.r[i+1] - self.r[i]

            phif = self.phi_flap(self.r[i])
            phie = self.phi_edge(self.r[i])

            Qf += FN_interp[i] * phif * dr
            Qe += FT_interp[i] * phie * dr

        return np.array([Qf, Qe])

    # ========================================================
    # NATURAL FREQUENCIES
    # ========================================================

    def natural_frequencies(self):

        if self.M is None:
            self.compute_M()

        if self.K is None:
            self.compute_K()

        omega_f = np.sqrt(self.K[0,0] / self.M[0,0])
        omega_e = np.sqrt(self.K[1,1] / self.M[1,1])

        f_f = omega_f / (2*np.pi)
        f_e = omega_e / (2*np.pi)

        return f_f, f_e


# ============================================================
# BUILD MODEL
# ============================================================

model = StructuralModel(
    r,
    m,
    EI_flap,
    EI_edge
)

# ============================================================
# COMPUTE MATRICES
# ============================================================

M = model.compute_M()
K = model.compute_K()
C = model.compute_C()

print("Mass matrix:")
print(M)

print("\nStiffness matrix:")
print(K)

print("\nDamping matrix:")
print(C)

# ============================================================
# NATURAL FREQUENCIES
# ============================================================

f_flap, f_edge = model.natural_frequencies()

print("\nFirst flapwise frequency [Hz]:", f_flap)
print("First edgewise frequency [Hz]:", f_edge)

# ============================================================
# GENERALIZED FORCES
# ============================================================

Q = model.generalized_forces(Rx, FN, FT)

print("\nGeneralized force vector:")
print(Q)

# ============================================================
# PLOTS
# ============================================================

plt.figure(figsize=(10,5))

plt.plot(Rx, FN, 'r-o', label='FN')
plt.plot(Rx, FT, 'b-o', label='FT')

plt.xlabel('Radius [m]')
plt.ylabel('Load [N/m]')
plt.grid(True)
plt.legend()

plt.title('BEM Loads')

plt.show()
