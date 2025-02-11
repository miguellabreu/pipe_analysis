from ansys.mapdl.core import launch_mapdl
import numpy as np

# Parameters
rho_s = 1000.0  # riser density 1000 kg/m3
L_s = 0.200     # riser length 200mm
De_s = 0.010    # riser outer diameter 10mm
Thic_riser = De_s / 2.0
Di_s = 0.0      # riser inner diameter
Dhydro = De_s   # riser hydrodynamic diameter
E_s = 0.5e9     # riser Young modulus
T_s = 0.0       # riser tension

# Dimensionless Variables (9)
ep = 0.02       # fluid damping in line-flow wake
eq = 0.04       # fluid damping in cross-flow wake

Ap = 96.0       # coupling coefficient line-flow wake
Aq = 12.0       # coupling coefficient cross-flow wake

Ca = 1.0        # added mass coefficient, fixed
Cm = Ca + 1.0
C0D = 1.2       # drag mean coefficient
C0L = 0.30      # lift mean coefficient
CiD0 = 0.2      # drag coefficient by fluctuating induced by the wake
CiL0 = 0.00
CT = 0.0

# Fluid Variables
rho_f = 1000.0  # 1000 kg/m3
phi_f = 0.0     # angle between the fluid velocity and the cross section in degrees
U_f = 0.05      # fluid velocity 5cm/s
pi = 3.1416

Ir = (pi * De_s**4) / 64.0 - (pi * Di_s**4) / 64.0  # Riser Inertia (radial plane)
Depth_sea = 1.0
A_s = 0.25 * pi * (De_s**2 - Di_s**2)

grav_z = 0.0

St = 0.2        # Strouhal number
freq_f = St * U_f / Dhydro
Tp_f = 1.0 / freq_f
omega_f = 2.0 * pi * freq_f * np.cos(phi_f)

# Begin Modeling
mapdl = launch_mapdl()

mapdl.prep7()

matpipe = 1
mapdl.mp('DENS', matpipe, rho_s)
mapdl.mp('EX', matpipe, E_s)
mapdl.mp('PRXY', matpipe, 0.3)

matwat = 2  # material number id of the ocean
mapdl.mp('DENS', matwat, rho_f)

mapdl.et(1, 'PIPE288')
mapdl.sectype(1, 'PIPE', 'riser')
mapdl.secdata(De_s, Thic_riser, 12, 0, 1, 0, 0, 0)
mapdl.keyopt(1, 4, 2)  # Thick pipe

mapdl.k(1, 0, 0, -L_s)
mapdl.k(2, 0, 0, 0)
mapdl.l(1, 2)

Esize_riser = L_s

mapdl.lesize(1, Esize_riser, 1)
mapdl.lmesh(1)

# Solution
mapdl.slashsolu()

mapdl.acel(0, 0, grav_z)  # gravity

mapdl.d(1, 'ALL')  # Plet BC's

mapdl.antype(4)  # Transient Analysis
mapdl.nlgeom(0)  # Large displacement

mapdl.trnopt('FULL')
mapdl.lumpm(1)  # Lumped mass
mapdl.trnopt('FULL', '', '', '', '', 'NMK')

time_tot = 1.4 * Tp_f
dt = 0.001
nsteps = int(time_tot / dt)

damp1 = 0.0
damp2 = 0.0

newmark1 = 0.25
newmark2 = 0.50
beta = newmark1
gamma = newmark2

nimp = 100

p = 2.0  # initial condition (p) Drag van der Pol oscillator
q = 2.0  # initial condition (q) Lift van der Pol oscillator
dp = 0.0  # initial condition (dp/dt) Drag van der Pol oscillator
dq = 0.0  # initial condition (dq/dt) Drag van der Pol oscillator

Cp = 2.0 * ep * omega_f * (p**2 - 1.0)
Kp = 4.0 * omega_f**2
Cq = eq * omega_f * (q**2 - 1.0)
Kq = omega_f**2

Cdrag = C0D
Clift = 0.0  # C0L

velX = 0.0
velY = 0.0
acelX = 0.0
acelY = 0.0

U2rel = (U_f - 0.5 * velX)**2 + 0.5 * velY**2

Fp = Ap * (0.5 * acelX) / Dhydro  # External Force condition (Fp)
Fq = Aq * (0.5 * acelY) / Dhydro  # External Force condition (Fq)

dp2 = (Fp - Cp * dp - Kp * p)
dq2 = (Fq - Cq * dq - Kq * q)

Pvector = np.zeros((nsteps, 1))
Qvector = np.zeros((nsteps, 1))
dispX = np.zeros((nsteps, 1))
dispY = np.zeros((nsteps, 1))

time_it = 0.0

# Add a breakpoint here
mapdl.run('/eof')

# Open the GUI to see the model running
mapdl.open_gui()

# Begin Loop Time
for it in range(1, nsteps + 1):
    time_it += dt

    mapdl.deltim(dt, 0, 0)
    mapdl.outres('ERASE')
    mapdl.outres('ALL', -1)
    mapdl.autots(0)
    mapdl.alphad(damp1)
    mapdl.betad(damp2)
    mapdl.time(time_it)
    mapdl.tintp('', newmark1, newmark2)
    mapdl.tintp(0.0)
    mapdl.kbc(1)  # stepped

    # Begin.VIV Drag oscillator
    Cp = 2.0 * ep * omega_f * (p**2 - 1.0)
    Fp = Ap * (0.5 * acelX) / Dhydro

    p_pred = p + dt * dp + 0.5 * dt * dt * (1.0 - 2.0 * beta) * dp2  # p predictor
    dp_pred = dp + 0.5 * dt * dt * (1.0 - 2.0 * beta) * dp2  # dp/dt predictor

    dp2 = (Fp - Cp * dp_pred - Kp * p_pred)  # d2p/dt2 (n+1)

    p = p_pred + beta * dt * dt * dp2
    dp = dp_pred + gamma * dt * dp2
    Pvector[it - 1, 0] = p

    # End.VIV Drag oscillator

    # Begin.VIV Lift oscillator
    Cq = eq * omega_f * (q**2 - 1.0)
    Kq = omega_f**2
    Fq = Aq * (0.5 * acelY) / Dhydro

    q_pred = q + dt * dq + 0.5 * dt * dt * (1.0 - 2.0 * beta) * dq2  # p predictor
    dq_pred = dq + 0.5 * dt * dt * (1.0 - 2.0 * beta) * dq2  # dp/dt predictor

    dq2 = (Fq - Cq * dq_pred - Kq * q_pred)  # d2p/dt2 (n+1)

    q = q_pred + beta * dt * dt * dq2
    dq = dq_pred + gamma * dt * dq2
    Qvector[it - 1, 0] = q
    # End.VIV Lift oscillator

    mapdl.solve()

    velX = mapdl.get_value('NODE', 2, 'VX')
    velY = mapdl.get_value('NODE', 2, 'VY')

    U2rel = (U_f - 0.5 * velX)**2 + 0.5 * velY**2
    Flift = -0.5 * 0.5 * q * C0L * rho_f * Dhydro * U2rel
    Fdrag = 0.5 * (C0D + 0.5 * p * CiD0) * rho_f * Dhydro * U2rel

    mapdl.sfbeam(1, 4, 'PRES', Fdrag, Fdrag)
    mapdl.sfbeam(1, 5, 'PRES', Flift, Flift)

    acelX = mapdl.get_value('NODE', 2, 'AX')
    acelY = mapdl.get_value('NODE', 2, 'AY')

    dispX[it - 1, 0] = mapdl.get_value('NODE', 2, 'UX') / De_s
    dispY[it - 1, 0] = mapdl.get_value('NODE', 2, 'UY') / De_s

# End Loop Time
mapdl.save()

# End Modeling