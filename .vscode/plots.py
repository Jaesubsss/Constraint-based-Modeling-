import numpy as np
import matplotlib.pyplot as plt

# Define x1 and x2 ranges
x1_values = np.linspace(0, 500, 100)
x2_values = np.linspace(0, 500, 100)

# Define the inequalities
inequality1 = (1500 - 4*x1_values) / 5
inequality2 = (1575 - 5*x1_values) / 3
inequality3 = (420 - x1_values) / 2

# Plot the inequalities
plt.plot(x1_values, inequality1, label='4x1 + 5x2 = 1500')
plt.plot(x1_values, inequality2, label='5x1 + 3x2 = 1575')
plt.plot(x1_values, inequality3, label='x1 + 2x2 = 420')

# Fill the feasible region
feasible_region = np.minimum(inequality1, np.minimum(inequality2, inequality3))
plt.fill_between(x1_values, 0, feasible_region, where=(feasible_region > 0), color='gray', alpha=0.5)

# Set labels and title
plt.xlabel('x1')
plt.ylabel('x2')
plt.title('Graphical Representation of Constraints')
plt.grid(True)
plt.legend()

# Show the plot
plt.show()