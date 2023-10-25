from event import Event
import sys
import numpy as np
from array import array
from ROOT import TFile, TTree

class Writer:

    def __init__(self, event):
        self.event = event

        self.fout = TFile('output.root', 'RECREATE')
        self.initOutputTree()

    def initOutputTree(self):
        self.tout = TTree('Sim', 'Sim') # output tree
        self.E_depoTotal = array('f', [0])
        # self.E_depoTotal = np.zeros((1,), dtype=np.float32)
        self.tout.Branch('E_depoTotal', self.E_depoTotal, 'E_depoTotal/F')

    def Write(self):
        self.fout.cd()

        for i in range(self.event.nEntry):
            self.event.Jump(i)

            self.E_depoTotal[0] = np.sum([depo.GetEnergyDeposit() for depo in self.event.depos]) 

            self.tout.Fill()

        self.tout.Write()

if __name__ == "__main__":
    event = Event(sys.argv[1])  
    w = Writer(event)
    w.Write()