from event import Event
from plotter import Plotter
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Look at first 10 events
n = 10
counter = 0
event = Event("/pnfs/dune/persistent/users/weishi/FD3/LArBath/Marley_Edepsim_noSecondaryDeposit/nue/edep_nue_35MeV_1kevts.root", 'Marley')
p = Plotter(event)

print("Look at first", n, "events")

#dx_all_evts = np.array([])
#trklength_all_evts = np.array([])
p.Jump(0, 0)

while counter < n:
    #if counter == 6 or counter == 5 or counter == 2 or counter == 36 or counter == 51:
    #if counter == 51:
    p.Draw('yz', value='time', energy= 'MeV')
    p.Draw('yz', value='charge', energy= 'MeV')
    #p.Draw('xz', value='time', energy= 'MeV')
    p.hist_LowE_dEdx()
    #p.hist_dx()
    #p.hist_trklength()

    # all events single edep length
    #dx_all_evts = np.concatenate((dx_all_evts, p.hist_dx()))
    # all events tracks length
    #trklength_all_evts = np.concatenate((trklength_all_evts, p.hist_trklength()))

    event.PrintVertex()
    event.PrintTracks(0, 5)
    p.Next(0) # dQ threshold

    counter += 1

"""plt.hist(dx_all_evts, range=(0,4), bins=100)
plt.xlabel('dx [cm]')
plt.draw()
plt.savefig('plots/single_edep_dx_all_evts.pdf')
plt.clf() # important to clear figure
plt.close()

plt.hist(trklength_all_evts, range=(0,6), bins=100)
plt.xlabel('track length [cm]')
plt.yscale("log")
plt.draw()
plt.savefig('plots/primary_track_length_all_evts.pdf')
plt.clf() # important to clear figure
plt.close()"""
