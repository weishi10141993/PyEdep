
## Description
Python class to Read/Write/Plot the [edep-sim](https://github.com/ClarkMcGrew/edep-sim/tree/master/io) files.

## Example Jupyter Notebook
```py
from event import Event
from plotter import Plotter
%matplotlib inline

event = Event("path/to/rootfile")
p = Plotter(event)

p.Draw('yz')

p.Next()
p.Draw('yz')

event.PrintVertex()
event.PrintTracks(0,8) 
```