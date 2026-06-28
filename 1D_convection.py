'''
CFD @ FIE
A.Y. 2025-2026
Professor: Edgardo A. Serafin
Student: Andrea Grimaldi

Final exam
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
from time import time

t0 = time()

#########################################################################################################################
#########################################################################################################################
# Starting parameters
t1 = 0.3
t2 = 0.9
t3 = 1.5
t4 = 2.0
tfinal = 3.4
nvol = 5000
xA = 1.2

# Derived parameters
config = 1   # 1 for case 1, 2 for case 2
make_animation = 0 # 1 to generate the animated GIF, 0 to not generate it

if config == 1:
    L = 2
    dx = (2 * L)/nvol
    uL = 0
    uR = 0
else:
    L = 3
    dx = (2 * L)/nvol
    uL = 1
    uR = 1

assert abs(xA) <= L, "xA should be in the study domain"

xp = np.linspace(-L + dx/2, L - dx/2, nvol)
u = np.zeros(nvol + 2)

# Given condition
if config == 1:
    u[0] = uL
    u[nvol + 1] = uR
    for i in range(1, nvol + 1):
        u[i] = math.exp(-2 * (xp[i - 1])*(xp[i - 1]))
else:
    u[0] = uL
    u[nvol + 1] = uR
    for i in range(1, nvol + 1):
        if xp[i - 1] >= -2.9 and xp[i - 1] <= -2:
            u[i] = 2
        else:
            u[i] = 1

# Allocating space for final data presentation
u1 = np.zeros(nvol + 2)
u2 = np.zeros(nvol + 2)
u3 = np.zeros(nvol + 2)
u4 = np.zeros(nvol + 2)
u_max = []
u_a = []
iA = np.argmin(np.abs(xp - xA)) + 1
times = []
snap = []
snap_freq = 20


#########################################################################################################################
#########################################################################################################################
# Helper functions
def getTimeStep(uMax):
    C = 0.8
    assert C <= 1, "CFL must be less or equal than 1"
    return (C * dx)/max(uMax, 1e-12)

def timeToNextTarget(t):
    if t < t1:
        return t1 - t
    if t < t2:
        return t2 - t
    if t < t3:
        return t3 - t
    if t < t4:
        return t4 - t
    return tfinal - t


#########################################################################################################################
#########################################################################################################################
# Loop
t = 0.0
step = 0

while t < tfinal:
    dt = getTimeStep(np.max(np.abs(u[1:-1])))
    dt = min(dt, timeToNextTarget(t))
    times.append(t)
    u_max.append(np.max(u[1:-1]))
    u_a.append(u[iA])

    ue = (u[1:-1] + u[2:])/2
    uw = (u[:-2] + u[1:-1])/2
    F = (u ** 2)/2
    Fe = np.where(ue > 0, F[1:-1], F[2:])
    Fw = np.where(uw > 0, F[:-2], F[1:-1])

    u[1:-1] = u[1:-1] - dt / dx * (Fe - Fw)

    t += dt

    if step % snap_freq == 0:
        snap.append((t, u[1:-1].copy()))

    step += 1

    if abs(t - t1) < 1e-9: u1[:] = u
    if abs(t - t2) < 1e-9: u2[:] = u
    if abs(t - t3) < 1e-9: u3[:] = u
    if abs(t - t4) < 1e-9: u4[:] = u


#########################################################################################################################
#########################################################################################################################
# u(x) over time
plt.figure()
plt.plot(xp, u1[1:-1], label=f"t = {t1}")
plt.plot(xp, u2[1:-1], label=f"t = {t2}")
plt.plot(xp, u3[1:-1], label=f"t = {t3}")
plt.plot(xp, u4[1:-1], label=f"t = {t4}")
plt.xlabel("x")
plt.ylabel("u")
plt.title(f"Burgers — Case {config}: u(x) at selected times")
plt.legend()
plt.grid(True)
plt.show()

# max(u(x)) over time
plt.figure()
plt.plot(times, u_max)
plt.xlabel("t")
plt.ylabel("max u")
plt.title(f"Burgers — Case {config}: maximum of u over time")
plt.grid(True)
plt.show()

# u(xA) over time
plt.figure()
plt.plot(times, u_a)
plt.xlabel("t")
plt.ylabel("u(xA)")
plt.title(f"Burgers — Case {config}: u({xA}) over time")
plt.grid(True)
plt.show()

if make_animation == 1:
    # u(x) for the study domain over time
    fig, ax = plt.subplots()
    line, = ax.plot([], [], lw=1.5)
    ax.set_xlim(-L, L)
    ax.set_ylim((0.9 if config == 2 else -0.1), (2.1 if config == 2 else 1.1))
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.grid(True)
    title = ax.set_title("")


    def init():
        line.set_data([], [])
        return line, title


    def animate(frame):
        t_f, u_f = frame
        line.set_data(xp, u_f)
        title.set_text(f"Burgers — Case {config}:  t = {t_f:.3f}")
        return line, title


    anim = animation.FuncAnimation(fig, animate, frames=snap,
                                   init_func=init, interval=40, blit=False)
    anim.save(f"burgers_case{config}.gif", writer="pillow", fps=25)
    print("animation saved")


#########################################################################################################################
#########################################################################################################################
t_end = time()
tCPU = t_end - t0
print('CPU time = {:1.2f} sec for {:1.0f} volumes'.format(tCPU, nvol))



