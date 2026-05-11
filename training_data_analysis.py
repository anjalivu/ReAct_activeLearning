import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data_0900 = pd.read_csv('training_data_0-90-0.csv')

t_FRP = data_0900['t_FRP'].values
w_f = data_0900['w_f'].values
prestress = data_0900['prestress'].values
E1_FRP = data_0900['E1_FRP'].values
pull_corners_needed = data_0900['pull_corners_needed'].values
has_inflection = data_0900['has_inflection'].values
radius = data_0900['radius'].values

EI = E1_FRP * w_f * 3*t_FRP**3 / 12
J = E1_FRP * w_f * t_FRP * (t_FRP**2 + w_f**2)

valid_mask = ~np.isnan(radius)
empty_mask = np.isnan(radius)



# Create 4 subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# t_FRP vs radius
axes[0, 0].scatter(t_FRP, radius, alpha=0.6)
axes[0, 0].set_xlabel('t_FRP')
axes[0, 0].set_ylabel('radius')
axes[0, 0].set_title('t_FRP vs radius')
axes[0, 0].grid(True, alpha=0.3)

# w_f vs radius
axes[0, 1].scatter(w_f, radius, alpha=0.6)
axes[0, 1].set_xlabel('w_f')
axes[0, 1].set_ylabel('radius')
axes[0, 1].set_title('w_f vs radius')
axes[0, 1].grid(True, alpha=0.3)

# prestress vs radius
axes[1, 0].scatter(prestress, radius, alpha=0.6)
axes[1, 0].set_xlabel('prestress')
axes[1, 0].set_ylabel('radius')
axes[1, 0].set_title('prestress vs radius')
axes[1, 0].grid(True, alpha=0.3)

# E1_FRP vs radius
axes[1, 1].scatter(E1_FRP, radius, alpha=0.6)
axes[1, 1].set_xlabel('E1_FRP')
axes[1, 1].set_ylabel('radius')
axes[1, 1].set_title('E1_FRP vs radius')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show



fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 5))

# First subplot: EI vs prestress
scatter1 = ax1.scatter(EI[valid_mask], prestress[valid_mask], 
                       c=radius[valid_mask], cmap='jet', s=500, alpha=0.9)
scatter_empty1 = ax1.scatter(EI[empty_mask], prestress[empty_mask], 
                             c='r', marker='x', s=100, label='Not Tristable')
fig.colorbar(scatter1, ax=ax1, label='radius')
ax1.set_xlabel('EI')
ax1.set_ylabel('prestress')
ax1.set_title('Bending stiffness vs prestress (colored by radius)')
ax1.grid(True, alpha=0.3)
ax1.legend()

# Second subplot: J vs prestress
scatter2 = ax2.scatter(J[valid_mask], prestress[valid_mask], 
                       c=radius[valid_mask], cmap='jet', s=500, alpha=0.9)
scatter_empty2 = ax2.scatter(J[empty_mask], prestress[empty_mask], 
                             c='r', marker='x', s=100, label='Not Tristable')
fig.colorbar(scatter2, ax=ax2, label='radius')
ax2.set_xlabel('J')
ax2.set_ylabel('prestress')
ax2.set_title('Torsional stiffness vs prestress (colored by radius)')
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.show()