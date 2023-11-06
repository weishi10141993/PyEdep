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

        self.vtx_xs = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('vtx_xs', self.vtx_xs, 'vtx_xs/F')

        self.vtx_proc = array('i', [0]) # genie process
        self.T_out.Branch('vtx_proc', self.vtx_proc, 'vtx_proc/I')

        self.vtx_nucl = array('i', [0]) # genie process
        self.T_out.Branch('vtx_nucl', self.vtx_nucl, 'vtx_nucl/I')

        self.E_nu = array('f', [0]) # energy availabe from vertex interaction (excluding energy lost inside nuclei)
        self.T_out.Branch('E_nu', self.E_nu, 'E_nu/F')

        self.E_avail = array('f', [0]) # energy availabe from vertex interaction (excluding energy lost inside nuclei)
        self.T_out.Branch('E_avail', self.E_avail, 'E_avail/F')    

        self.E_list = np.zeros((6,), dtype=np.float32) # E avail for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_list', self.E_list, 'E_list[6]/F') 

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

            # proc = str(self.event.info['vtx_proc']) + '-' + str(self.event.info['vtx_nucl'])
            # v = self.stat.setdefault(proc, 0)
            # self.stat[proc] = v + 1

            self.vtx_xs[0] = self.event.info['vtx_xs']
            self.vtx_proc[0] = self.event.info['vtx_proc']
            self.vtx_nucl[0] = self.event.info['vtx_nucl']

            self.E_nu[0] = self.GetEnu()
            self.E_depoTotal[0] = self.event.info['E_depoTotal']
            self.E_avail[0] = self.event.info['E_avail']
            self.E_list[:] = self.event.info['E_list']
            self.E_depoList[:] = self.event.info['E_depoList']

            self.T_out.Fill()

        self.T_out.Write()
        # print(self.stat)
    
    def GetEnu(self):
        filename = self.event.GetFileName().split('/')[-1]
        filename = filename.split('_')[-2]
        filename = filename.replace('GeV', '')
        return float(filename)

if __name__ == "__main__":
    if (len(sys.argv)>2):
        outfile = sys.argv[2]
    else:
        outfile = 'output.root'
    event = Event(sys.argv[1])  
    w = Writer(event, outfile)
    w.Write()