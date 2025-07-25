import matplotlib.pyplot as plt

"""
Draw a plot for complexity of the dual approach as a function of allowed oracle calls
"""
""" OLD
"""
# k = 76, depth 2
C = [152.64, 150.26, 146.28, 141.75, 140.99, 138.89, 135.36, 130.94]

M = [14.68, 15.88, 16.48, 18.32, 19.9, 20.82, 23.24, 28.32]

PARAMS = [[19, 18, 10], [19, 18, 9], [17, 17, 8], [17, 17, 7], [16, 15, 6], [17, 17, 6], [16, 15, 5], [16, 15, 4]]
# Plot
plt.figure(figsize=(8, 6))
plt.scatter(M, C, color='b', marker='o')  # Scatter plot of M vs. C


# Adding labels from PARAMS array
for i, (delta1, delta2, t) in enumerate(PARAMS):
    label = f"({delta1}, {delta2}, {t})"
    plt.text(M[i], C[i], label, fontsize=9, ha='left')  # Customize label placement with ha='right'
# Adding horizontal line at y = 143
plt.axhline(y=143, color='red', linestyle='-', linewidth=0.5, label = 'shifted bjmm')  # Thin red line at y = 143

# Adding axis labels and title
plt.legend()
plt.xlabel('$\log 2$ (Oracle calls)')
plt.ylabel('complexity')
plt.title('$k = 76$, $p = 127$, $E = \{ 2^i, i = 0,..., 6 \} $ ')
plt.savefig('k76_depth2.png', format='png', dpi=300)
# Display the plot
#plt.show()


# k = 76, depth 3
C = [151.09, 146.59, 140.99, 139.85, 135.44, 134.48, 128.74]
M = [19.67, 20.74, 21.49, 22.77, 23.88, 26.14, 27.49]
PARAMS = [[15, 14, 14, 6], [16, 16, 14, 6], [17, 16, 16, 6], [15, 14, 14, 5], [16, 15, 15, 5], [14, 13, 13, 4], [15, 14, 14, 4]]

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(M, C, color='b', marker='o')  # Scatter plot of M vs. C


# Adding labels from PARAMS array
for i, (delta1, delta2, delta3, t) in enumerate(PARAMS):
    label = f"({delta1}, {delta2}, {delta3}, {t})"
    plt.text(M[i], C[i], label, fontsize=9, ha='left')  # Customize label placement with ha='right'
# Adding horizontal line at y = 143
plt.axhline(y=143, color='red', linestyle='-', linewidth=0.5, label = 'shifted bjmm')  # Thin red line at y = 143

# Adding axis labels and title
plt.legend()
plt.xlabel('$\log 2$ (Oracle calls)')
plt.ylabel('complexity')
plt.title('$k = 76$, $p = 127$, $E = \{ 2^i, i = 0,..., 6 \} $ ')
plt.savefig('k76_depth3.png', format='png', dpi=300)
# Display the plot
#plt.show()

# k = 50, depth 3, bigger E
C = [143.52, 135.33, 127.89, 119.78, 112.3, 105.04]
M = [17.78, 18.45, 20.93, 22.21, 26.87, 28.57]
PARAMS = [[12, 10, 10, 5], [12, 12, 11, 5], [11, 11, 10, 4], [12, 12, 11, 4], [11, 11, 10, 3], [12, 12, 11, 3]]

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(M, C, color='b', marker='o')  # Scatter plot of M vs. C


# Adding labels from PARAMS array
for i, (delta1, delta2, delta3, t) in enumerate(PARAMS):
    label = f"({delta1}, {delta2}, {delta3}, {t})"
    plt.text(M[i], C[i], label, fontsize=9, ha='left')  # Customize label placement with ha='right'
# Adding horizontal line at y = 143
plt.axhline(y=143, color='red', linestyle='-', linewidth=0.5, label = 'shifted bjmm')  # Thin red line at y = 143

# Adding axis labels and title
plt.legend()
plt.xlabel('$\log 2$ (Oracle calls)')
plt.ylabel('complexity')
plt.title('$k = 50$, $p = 127$, $E = \{ (-2)^i, i = 0,..., 13 \} $ ')
plt.savefig('k76_depth3_z=14.png', format='png', dpi=300)
# Display the plot
#plt.show()
