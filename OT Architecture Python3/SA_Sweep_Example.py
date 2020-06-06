#!/usr/bin/env python

from opentire import OpenTire
from opentire.Core import TireState
from opentire.Core import TIRFile

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    # Initialize the tire model
    openTire = OpenTire()
    myTireModel = openTire.createmodel('PAC2002')

    # Initialize the tire state
    state = TireState()
    state['FZ'] = 1500
    state['IA'] = 0.0
    state['SR'] = 0.0
    state['SA'] = 0.0
    state['FY'] = 0.0
    state['V'] = 10.0
    state['P'] = 260000

    # Define the slip angle range
    slip_angles = np.linspace(-15, 15, 50) * 3.14 / 180

    # Print out some pretty formatting
    print('OpenTire Slip Angle Sweep Demo\n')
    print('{0:>10} | {1:>10} | {2:>10} | {3:>10} | {4:>10}'
          .format('SA [deg]',
                  'FZ [N]',
                  'FY [N]',
                  'MZ [Nm]',
                  'MX [Nm]'))
    print('=' * 62)

    # Pre-allocate an array for the result
    FyArray = []
    MzArray = []
    MxArray = []

    # Calculate and print out the tire model outputs
    for sa in slip_angles:
        state['SA'] = sa
        myTireModel.solve(state)

        # Grab results for plotting
        FyArray.append(state['FY'])
        MzArray.append(state['MZ'])
        MxArray.append(state['MX'])

        print('{0:>10.1f} | {1:>10.0f} | {2:>10.1f} | {3:>10.1f} | {4:>10.1f}'
              .format(state['SA'] * 180 / 3.14,
                      state['FZ'],
                      state['FY'],
                      state['MZ'],
                      state['MX']))

# Plot results
fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.plot(slip_angles, FyArray)
ax1.set_xlabel("Slip Angle (rad)")
ax1.set_ylabel("FY (N)")
ax1.set_title("Lateral Force")

ax2 = fig.add_subplot(223)
ax2.plot(slip_angles, MzArray)
ax2.set_xlabel("Slip Angle (rad)")
ax2.set_ylabel("MZ (Nm)")
ax2.set_title("Self aligning moment")

ax3 = fig.add_subplot(224)
ax3.plot(slip_angles, MxArray)
ax3.set_xlabel("Slip Angle (rad)")
ax3.set_ylabel("MX (Nm)")
ax3.set_title("Overturning  moment")
plt.tight_layout()
plt.show()
