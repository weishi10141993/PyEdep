
## Description
Python class to Read/Write/Plot the [edep-sim](https://github.com/ClarkMcGrew/edep-sim/tree/master/io) files.

## Example Jupyter Notebook
```py
from event import Event
from plotter import Plotter
%matplotlib inline

event = Event("path/to/rootfile")
p = Plotter(event)

p.Draw('yz', value='time')
p.Draw('yz', value='charge')

p.Next()
p.Draw('yz')

event.PrintVertex()
event.PrintTracks(0, 8)
event.PrintTrack(0) 
```

## Writer
Summarize list of edep-sim root files into a customized output tree
```py
python3 writer.py 'input_file*.root' output_file.root
```

The branches are defined as follows:
| Branch        | unit           | description  |
| :------------ |:-------------| :-----|
| nu_pdg        | int | neutrino pdg code |
| nu_xs      | cm^2      | interaction cross section |
| nu_proc | int | 10+X: CC; 20+X: NC; X: 1: QES; 2: RES; 3: DIS; 4: COH, 5: MEC  |
| nu_nucl | int | interact with proton or neutron |
| E_nu | MeV | true neutrino energy |
| E_avail | MeV | available energy from initial state particles, including mass if meson or lepton |
| E_availList | array, MeV | 0: mu/e; 1: proton; 2: neutron; 3: pi+/-; 4: pi0; 5: others |
| E_depoTotal | MeV | total deposited energy  |
| E_depoList | array, MeV | similar to E_availList but for deposited energy including all children |
