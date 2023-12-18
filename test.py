from event import Event
from plotter import Plotter

# Look at first 10 events
n = 10
counter = 0
event = Event("edep_nue_5MeV_1kevts.root", 'Marley')
p = Plotter(event)

print("Look at first", n, "events")

while counter < n:
  print("event", counter)

  p.Draw('yz', value='time')
  p.Draw('yz', value='charge')
  p.hist_dqdx()
  event.PrintVertex()
  event.PrintTracks(0, 5)
  p.Next()

  counter += 1
