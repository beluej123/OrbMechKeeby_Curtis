import numpy as np
import rv_from_r0v0_ta

def rv_from_r0v0_ta_test():
    '''
    This program computes the state vector [R, V] from the initial
    state vector [R0, V0] and the change in true anomaly, using the data.

    mu - gravitational parameter (km^3/s^2)
    R0 - the initial position vector (km)
    V0 - the initial velocity vector (km/s)
    r0 - magnitude of R0
    v0 - magnitude of V0
    R  - final position vector (km)
    V  - final velocity vector (km/s)
    r  - magnitude of R
    v  - magnitude of V
    dt - change in true anomaly (degrees)

    User py-function required: rv_from_r0v0_ta
    '''
    mu = 398600.4418

    #...Input data
    R0 = np.array([8182.4, -6865.9, 0])
    V0 = np.array([0.47572, 8.8116, 0])
    dt = 120
    #...End input data

    #...Algorithm 2.3:
    R, V = rv_from_r0v0_ta.rv_from_r0v0_ta(R0, V0, dt, mu)

    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)

    print("--------------------------------------------------------------------")
    print("\n Initial state vector:")
    print(f"\n r = [{R0[0]:.1f}, {R0[1]:.1f}, {R0[2]:.1f}] (km)")
    print(f" magnitude = {r0:.2f}")
    print(f"\n v = [{V0[0]:.5f}, {V0[1]:.4f}, {V0[2]:.1f}] (km/s)")
    print(f" magnitude = {v0:.4f}")
    print(f"\n\n State vector after {dt} degree change in true anomaly:")
    print(f"\n r = [{R[0]:.2f}, {R[1]:.2f}, {R[2]:.1f}] (km)")
    print(f" magnitude = {r:.2f}")
    print(f"\n v = [{V[0]:.5f}, {V[1]:.5f}, {V[2]:.1f}] (km/s)")
    print(f" magnitude = {v:.5f}")
    print("--------------------------------------------------------------------")

if __name__ == "__main__":
    rv_from_r0v0_ta_test()
