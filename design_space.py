import numpy as np

N = 50
theta_range = [-90, 90]
phi_range = [-90, 90]
t_frp_range = [0.025, 0.060]
w_f_range = [6.0, 15.0]
prestress_range = [4.0, 8.0]
E1_FRP_range = [50000.0, 250000.0]

n_dim = 6

seed = 69
np.random.seed(seed)

bins = np.arange(N)
samples = np.zeros((N, n_dim))
for i in range(n_dim):
    np.random.shuffle(bins)
    # (bin + random[0,1)) / N ensures values stay in [0, 1)
    samples[:, i] = (bins) / N
noise = np.random.uniform(0, 1/N, size=samples.shape)
samples = samples + noise


design_space = np.zeros((N, n_dim))

for i in range(N):
    sample = samples[i]
    theta = np.interp(sample[0], [0, 1], theta_range)
    phi = np.interp(sample[1], [0, 1], phi_range)
    t_FRP = np.interp(sample[2], [0, 1], t_frp_range)
    w_f = np.interp(sample[3], [0, 1], w_f_range)
    prestress = np.interp(sample[4], [0, 1], prestress_range)
    E1_FRP = np.interp(sample[5], [0, 1], E1_FRP_range)


    design_space[i] = [theta, phi, t_FRP, w_f, prestress, E1_FRP]

np.savetxt('design_space.csv', design_space, delimiter=',')


