from event import Event
import sys
import numpy as np
from array import array
from ROOT import TFile, TTree

class Writer:

    def __init__(self, event, outfile='output.root'):
        self.event = event

        self.f_out = TFile(outfile, 'RECREATE')
        self.initOutputTree()

    def initOutputTree(self):
        self.T_out = TTree('Sim', 'Sim') # output tree

        self.nu_pdg = array('i', [0]) # 
        self.T_out.Branch('nu_pdg', self.nu_pdg, 'nu_pdg/I')

        self.nu_xs = array('f', [0]) 
        self.T_out.Branch('nu_xs', self.nu_xs, 'nu_xs/F')

        self.nu_proc = array('i', [0]) # genie process
        self.T_out.Branch('nu_proc', self.nu_proc, 'nu_proc/I')

        self.nu_nucl = array('i', [0]) # which nucleon
        self.T_out.Branch('nu_nucl', self.nu_nucl, 'nu_nucl/I')

        self.E_nu = array('f', [0]) # true neutrino energy
        self.T_out.Branch('E_nu', self.E_nu, 'E_nu/F')

        self.E_avail = array('f', [0]) # energy availabe from vertex interaction (excluding energy lost inside nuclei)
        self.T_out.Branch('E_avail', self.E_avail, 'E_avail/F')    

        self.E_availList = np.zeros((6,), dtype=np.float32) # E avail for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_availList', self.E_availList, 'E_availList[6]/F') 

        self.E_depoTotal = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('E_depoTotal', self.E_depoTotal, 'E_depoTotal/F')

        self.E_depoList = np.zeros((6,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_depoList', self.E_depoList, 'E_depoList[6]/F')



    def Write(self):
        # self.stat = {
        # }
        self.f_out.cd()

        for i in range(self.event.nEntry):
        # for i in range(100):
            self.event.Jump(i)

            # proc = str(self.event.info['nu_proc']) + '-' + str(self.event.info['nu_nucl'])
            # v = self.stat.setdefault(proc, 0)
            # self.stat[proc] = v + 1
            self.nu_pdg[0] = self.event.info['nu_pdg']
            self.nu_xs[0] = self.event.info['nu_xs']
            self.nu_proc[0] = self.event.info['nu_proc']
            self.nu_nucl[0] = self.event.info['nu_nucl']
            self.E_nu[0] = self.event.info['E_nu']

            # self.E_nu[0] = self.GetEnu()
            self.E_depoTotal[0] = self.event.info['E_depoTotal']
            self.E_avail[0] = self.event.info['E_avail']
            self.E_availList[:] = self.event.info['E_availList']
            self.E_depoList[:] = self.event.info['E_depoList']

            self.T_out.Fill()

        self.T_out.Write()
        # print(self.stat)
    
    def GetEnu(self):
        # no need, can read from gRooTracker tree directly
        filename = self.event.GetFileName().split('/')[-1]
        filename = filename.split('_')[-2]
        filename = filename.replace('GeV', '')
        return float(filename)*1000  # to MeV

if __name__ == "__main__":
    if (len(sys.argv)>2):
        outfile = sys.argv[2]
    else:
        outfile = 'output.root'
    event = Event(sys.argv[1])  
    w = Writer(event, outfile)
    w.Write()