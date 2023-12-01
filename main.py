from assignement8 import *
from tqdm import tqdm
import numpy as np
from scipy.integrate import solve_ivp


# PLOT up to the first crossing for mu = 0.008
mu = 0.008
for i in range(len(vaps_L3(mu))):
    if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
        v0 = veps_L3(mu)[:, i].real
        print(v0)  # Eigvec for positive eigval
        fig = plt.axes()
        fig.set_aspect('equal', 'box')
        fig.scatter(L3(mu), 0, c='crimson')

        x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
        if x0[1] > 0:
            x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
        a = Poincar_Map(x0, 1, mu,1)
        l=[]
        data=np.linspace(0, a[1], 1000)
        for i in range(len(data)):
            l.append(data[i][0])
        sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
                        t_eval=np.array(l), rtol=3e-14, atol=1e-14)
        fig.plot(sol.y[0], sol.y[1], c='green')

plt.show()


###Plot for varying mu
# Intervals for the solutions
I1 = np.linspace(0.001, 0.015, round((0.015 - 0.001) / 0.00001) + 1)
I2 = np.linspace(0.015, 0.05, round((0.05 - 0.015) / 0.0001) + 1)
I3 = np.linspace(0.05, 0.49, round((0.49 - 0.05) / 0.001) + 1)
#
# # Function to compute all the crossings
def Calcul(mu):
    for i in range(len(vaps_L3(mu))):
        if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
            v0 = veps_L3(mu)[:, i].real  # Eigvec for positive eigval
    x0 = np.array([L3(mu), 0, 0, 0]) + 1E-6 * v0
    if x0[1] > 0:
        x0 = np.array([L3(mu), 0, 0, 0]) - 1E-6 * v0
    a = Poincar_Map(x0, 1, mu,1)[0]
    pbar.update(1)
    return np.array([a[0][0], a[1][0], a[2][0], a[3][0]])

## Computation of the crossings of the Poincaré section for the three intervals

# results_I1_x = np.array([])
# results_I1_y = np.array([])
# results_I1_xp = np.array([])
# results_I1_yp = np.array([])
# #
# tasks = I1
# pbar = tqdm(total=len(tasks))
# #
# for mu in tasks:
#     a = Calcul(mu)
#     results_I1_x = np.append(results_I1_x, a[0])
#     results_I1_y = np.append(results_I1_y, a[1])
#     results_I1_xp = np.append(results_I1_xp, a[2])
#     results_I1_yp = np.append(results_I1_yp, a[3])
# pbar.close()
#
# print(results_I1_x)
# np.savetxt('results_I1_x.csv', results_I1_x, delimiter=',')
# np.savetxt('results_I1_y.csv', results_I1_y, delimiter=',')
# np.savetxt('results_I1_xp.csv', results_I1_xp, delimiter=',')
# np.savetxt('results_I1_yp.csv', results_I1_yp, delimiter=',')
#
#
# results_I2_x = np.array([])
# results_I2_y = np.array([])
# results_I2_xp = np.array([])
# results_I2_yp = np.array([])
#
# tasks = I2
# pbar = tqdm(total=len(tasks))
#
# for mu in tasks:
#     a = Calcul(mu)
#     results_I2_x = np.append(results_I2_x, a[0])
#     results_I2_y = np.append(results_I2_y, a[1])
#     results_I2_xp = np.append(results_I2_xp, a[2])
#     results_I2_yp = np.append(results_I2_yp, a[3])
#
# pbar.close()
#
# np.savetxt('results_I2_x.csv', results_I2_x, delimiter=',')
# np.savetxt('results_I2_y.csv', results_I2_y, delimiter=',')
# np.savetxt('results_I2_xp.csv', results_I2_xp, delimiter=',')
# np.savetxt('results_I2_yp.csv', results_I2_yp, delimiter=',')

#
# results_I3_x = np.array([])
# results_I3_y = np.array([])
# results_I3_xp = np.array([])
# results_I3_yp = np.array([])
#
# tasks = I3
# pbar = tqdm(total=len(tasks))
#
# for mu in tasks:
#     a = Calcul(mu)
#     results_I3_x = np.append(results_I3_x, a[0])
#     results_I3_y = np.append(results_I3_y, a[1])
#     results_I3_xp = np.append(results_I3_xp, a[2])
#     results_I3_yp = np.append(results_I3_yp, a[3])
#
# pbar.close()
#
# np.savetxt('results_I3_x.csv', results_I3_x, delimiter=',')
# np.savetxt('results_I3_y.csv', results_I3_y, delimiter=',')
# np.savetxt('results_I3_xp.csv', results_I3_xp, delimiter=',')
# np.savetxt('results_I3_yp.csv', results_I3_yp, delimiter=',')

###First part: plot for the explanation of the firsts disconinuitiesE
# mu = 0.00925
# for i in range(len(vaps_L3(mu))):
#     if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
#         v0 = veps_L3(mu)[:, i].real
#         print(v0)  # Eigvec for positive eigval
#         # fig = plt.axes()
#         # fig.set_aspect('equal', 'box')
#         # fig.scatter(L3(mu), 0, c='crimson')
#
#         x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
#         if x0[1] > 0:
#             x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
#         a = Poincar_Map(x0, 1, mu,1)
#         l=[]
#         data=np.linspace(0, a[1], 1000)
#         for i in range(len(data)):
#             l.append(data[i][0])
#         sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
#                         t_eval=np.array(l), rtol=3e-14, atol=1e-14)
#         plt.plot(sol.y[0], sol.y[1], c='green',label="mu=0.00925")
#
# mu = 0.00926
#
# for i in range(len(vaps_L3(mu))):
#     if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
#         v0 = veps_L3(mu)[:, i].real
#         print(v0)  # Eigvec for positive eigval
#         # fig = plt.axes()
#         # fig.set_aspect('equal', 'box')
#         # fig.scatter(L3(mu), 0, c='crimson')
#
#         x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
#         if x0[1] > 0:
#             x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
#         a = Poincar_Map(x0, 1, mu,1)
#         l=[]
#         data=np.linspace(0, a[1], 1000)
#         for i in range(len(data)):
#             l.append(data[i][0])
#         sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
#                         t_eval=np.array(l), rtol=3e-14, atol=1e-14)
#         plt.plot(sol.y[0], sol.y[1], c='red',label="mu=0.00926")
#
# mu = 0.00924
# for i in range(len(vaps_L3(mu))):
#     if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
#         v0 = veps_L3(mu)[:, i].real
#         print(v0)  # Eigvec for positive eigval
#
#         x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
#         if x0[1] > 0:
#             x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
#         a = Poincar_Map(x0, 1, mu,1)
#         l=[]
#         data=np.linspace(0, a[1], 1000)
#         for i in range(len(data)):
#             l.append(data[i][0])
#         sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
#                         t_eval=np.array(l), rtol=3e-14, atol=1e-14)
#         plt.plot(sol.y[0], sol.y[1], c='pink',label="mu=0.00924")
#
# mu = 0.00926
# for i in range(len(vaps_L3(mu))):
#     if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
#         v0 = veps_L3(mu)[:, i].real
#         print(v0)  # Eigvec for positive eigval
#
#         x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
#         if x0[1] > 0:
#             x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
#         a = Poincar_Map(x0, 3, mu,1)
#         l=[]
#         data=np.linspace(0, a[1], 1000)
#         for i in range(len(data)):
#             l.append(data[i][0])
#         sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
#                         t_eval=np.array(l), rtol=3e-14, atol=1e-14)
#         s=len(sol.y[0])
#         plt.plot(sol.y[0], sol.y[1], linestyle='--',c='red')
# plt.xlim(0.6,1.1)
# plt.ylim(-0.10,0.10)
# plt.plot(sol.y[0],np.zeros(len(sol.y[0])))
# plt.legend()
# plt.show()

###Second part: plot for the explanation of the last discontinuity
# mu = 0.0198
# for i in range(len(vaps_L3(mu))):
#     if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
#         v0 = veps_L3(mu)[:, i].real
#         print(v0)  # Eigvec for positive eigval
#
#         x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
#         if x0[1] > 0:
#             x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
#         a = Poincar_Map(x0, 1, mu,1)
#         l=[]
#         data=np.linspace(0, a[1], 1000)
#         for i in range(len(data)):
#             l.append(data[i][0])
#         sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
#                         t_eval=np.array(l), rtol=3e-14, atol=1e-14)
#         s=len(sol.y[0])
#         plt.plot(sol.y[0], sol.y[1], linestyle='-',c='grey',label="mu=0.0198")
#
# mu = 0.0202
# for i in range(len(vaps_L3(mu))):
#     if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
#         v0 = veps_L3(mu)[:, i].real
#         print(v0)  # Eigvec for positive eigval
#
#         x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
#         if x0[1] > 0:
#             x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
#         a = Poincar_Map(x0, 1, mu,1)
#         l=[]
#         data=np.linspace(0, a[1], 1000)
#         for i in range(len(data)):
#             l.append(data[i][0])
#         sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
#                         t_eval=np.array(l), rtol=3e-14, atol=1e-14)
#         s=len(sol.y[0])
#         plt.plot(sol.y[0], sol.y[1], linestyle='-',c='green',label="mu=0.0202")
#
# mu = 0.0202
# for i in range(len(vaps_L3(mu))):
#     if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
#         v0 = veps_L3(mu)[:, i].real
#         print(v0)  # Eigvec for positive eigval
#
#         x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
#         if x0[1] > 0:
#             x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
#         a = Poincar_Map(x0, 2, mu,1)
#         l=[]
#         data=np.linspace(0, a[1], 1000)
#         for i in range(len(data)):
#             l.append(data[i][0])
#         sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
#                         t_eval=np.array(l), rtol=3e-14, atol=1e-14)
#         plt.plot(sol.y[0], sol.y[1], linestyle='-',c='green')
#
# mu = 0.0198
# for i in range(len(vaps_L3(mu))):
#     if np.isreal(vaps_L3(mu)[i]) and vaps_L3(mu)[i] > 0:
#         v0 = veps_L3(mu)[:, i].real
#         print(v0)  # Eigvec for positive eigval
#
#         x0 = np.array([L3(mu), 0, 0, 0]) + 10**-6 * v0
#         if x0[1] > 0:
#             x0 = np.array([L3(mu), 0, 0, 0]) - 10**-6 * v0
#         a = Poincar_Map(x0, 2, mu,1)
#         l=[]
#         data=np.linspace(0, a[1], 1000)
#         for i in range(len(data)):
#             l.append(data[i][0])
#         sol = solve_ivp(fun=lambda t, y: RTBP(t, y, mu), t_span=[0, a[1]], y0=x0,
#                         t_eval=np.array(l), rtol=3e-14, atol=1e-14)
#         plt.plot(sol.y[0], sol.y[1], linestyle='-',c='grey')
#
# plt.xlim(-1.20,-0.90)
# plt.ylim(-0.10,0.10)
# plt.plot(sol.y[0],np.zeros(len(sol.y[0])))
# plt.legend()
# plt.show()