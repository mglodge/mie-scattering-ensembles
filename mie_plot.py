'''

Simple code to read the output of "Mie Scattering (ensembles).c" and plot graphs of the results:

    - Phase function (total average scattering amplitude S_total vs angle)

    - Decomposed phase function (S1/S2/S_total vs angle)

    - Degree of linear polarisation (DOLP) vs angle)

                                Matt Lodge (17/07/23)

'''

import matplotlib.pyplot as plt
import pandas as pd

cols = ['angle', 'S1', 'S2', 'S_total', 'DOLP']
data = pd.read_csv('S1_and_S2_vs_angle.txt', names=cols)



# plot S total vs angle

plt.plot(data['angle'], data['S_total'], 'k', label='$|S1|^2+|S2|^2$ - Total')

plt.xlabel(r'$\theta$ (deg)', fontsize=13)
plt.ylabel('Intensity', fontsize=13)

plt.legend(prop={'size': 16})

# show graph
plt.yscale('log')
plt.show()



# plot S1, S2, and S total vs angle

plt.plot(data['angle'], data['S1'], 'b', label='$|S1|^2$ - Perpendicular')
plt.plot(data['angle'], data['S2'], 'g', label='$|S2|^2$ - Parallel')
plt.plot(data['angle'], data['S_total'], 'm', label='$|S1|^2+|S2|^2$ - Total')

plt.xlabel(r'$\theta$ (deg)', fontsize=13)
plt.ylabel('Intensity', fontsize=13)

plt.legend(prop={'size': 16})

# show graph
plt.yscale('log')
plt.show()


# plot degree of linear polarisation vs angle


plt.plot(data['angle'], data['DOLP'], 'b', label='D.O.L.P')

plt.xlabel(r'$\theta$ (deg)', fontsize=13)
plt.ylabel('% Degree of Linear Polarisation', fontsize=13)

plt.legend(prop={'size': 16})

# show graph
plt.yscale('linear')
plt.show()



