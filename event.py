import ROOT
from ROOT import TG4Event, TFile, TChain

import sys, os
import numpy as np

class Event:

    def __init__(self, fileName):
        self.fileName = fileName
        self.ReadTree()

        self.currentEntry = 0
        # self.Jump(self.currentEntry)

    # ------------------------
    def ReadTree(self):
        # self.rootFile = TFile(self.fileName)
        self.simTree = TChain("EDepSimEvents")
        self.simTree.Add(self.fileName)
        self.genieTree = TChain("DetSimPassThru/gRooTracker")
        self.genieTree.Add(self.fileName)
        # self.simTree = self.rootFile.Get("EDepSimEvents")
        # self.genieTree = self.rootFile.Get("DetSimPassThru/gRooTracker")

        self.nEntry = self.simTree.GetEntries()
        if self.nEntry != self.genieTree.GetEntries():
            print("Edep-sim tree and GENIE tree number of entries do not match!")
            sys.exit()        
        
        self.event = TG4Event()
        self.simTree.SetBranchAddress("Event", self.event)
    
    # ------------------------    
    def Jump(self, entryNo):
        print(f'reading event {entryNo}/{self.nEntry}')

        self.currentEntry = entryNo
        self.simTree.GetEntry(entryNo)

        self.ReadVertex()
        self.ReadTracks()
        self.ReadEnergyDepo('SimEnergyDeposit')

        self.info = {}
        self.info['E_nu'] = 0
        self.info['E_avail'] = 0
        self.info['E_availList'] = np.zeros(6) # lepton, proton, neutron, pi+-, pi0, others.
        self.info['E_depoTotal'] = 0
        self.info['E_depoList'] = np.zeros(6) # lepton, proton, neutron, pi+-, pi0, others.
        self.info['nu_pdg'] = 0
        self.info['nu_xs'] = self.vertex.GetCrossSection()
        self.info['nu_proc'], self.info['nu_nucl'] = self.GetReaction()
        self.FillEnergyInfo()
        self.ReadGenie()

    # ------------------------
    def ReadGenie(self):
        # gRooTracker info
        # https://github.com/GENIE-MC/Generator/blob/master/src/Apps/gNtpConv.cxx#L1837
        # StdHepStatus: 0: initial state; 1: final state particles; others: intermediate transport
        # following assumes the first particle is always the neutrino. 
        self.genieTree.GetEntry(self.currentEntry)
        self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
        self.info['E_nu'] = self.genieTree.StdHepP4[3]*1000

    # ------------------------
    def ReadVertex(self):
        primaries = np.array(self.event.Primaries)
        if (primaries.size != 1):
            print("Number of primaries not equal to 1 (not neutrino vertex)!")
            return

        self.vertex = primaries[0]

    #--------------------------
    def GetReaction(self):
        txt_list = self.vertex.GetReaction().split(';')
        proc = ''
        nucl = 0
        for x in txt_list[2:]:
            if 'proc:' in x:
                proc = x.replace('proc:', '').replace('Weak[','').replace('],', '')
            elif 'N:' in x:
                nucl = int(x.replace('N:', ''))
        # print(proc, nucl)
        proc_num = 0
        if (proc.startswith('CC')): 
            proc_num = 10
            proc = proc.replace('CC', '')
        elif (proc.startswith('NC')): 
            proc_num = 20
            proc = proc.replace('NC', '')
        else: 
            proc_num = 30
        proc_dict = {
            'QES' : 1,
            'RES' : 2,
            'DIS' : 3,
            'COH' : 4,
            'MEC' : 5,
        }
        proc_num += proc_dict.get(proc, 0)
        return proc_num, nucl

    # ------------------------
    def ReadTracks(self):
        self.tracks = np.array(self.event.Trajectories)
        for i in range(self.tracks.size):
            self.tracks[i].energy = {}
            self.tracks[i].association = {}
            self.tracks[i].energy['depoTotal'] = 0
            self.tracks[i].association['depoList'] = []
            self.tracks[i].association['children'] = []
            self.tracks[i].association['ancestor'] = i
            
        for i in range(self.tracks.size):
            track = self.tracks[i]
            parId = track.GetParentId()
            if parId == -1: continue
            self.tracks[parId].association['children'].append(i)
            while parId != -1:
                track.association['ancestor'] = parId
                parId = self.tracks[parId].GetParentId()

    # ------------------------
    def ReadEnergyDepo(self, detName):
        self.depos = np.array(self.event.SegmentDetectors[detName])

        # add depo info to tracks
        # depoList = self.FindDepoListFromTrack(1)
        # depoEnergy = np.sum([depo.GetEnergyDeposit() for depo in self.depos[depoList]])
        # print('debug: muon deposit energy:', depoEnergy)

        for i, depo in enumerate(self.depos):      
            trkId = depo.Contrib[0]
            edep = depo.GetEnergyDeposit()
            # edepSecond = depo.GetSecondaryDeposit()
            # trkLength = depo.GetTrackLength()
            track = self.tracks[trkId]
            track.association['depoList'].append(i)
            track.energy['depoTotal'] += edep

        # E_tot = np.sum([depo.GetEnergyDeposit() for depo in self.depos])
        # print('total deposit energy: ', E_tot)

    # ------------------------
    def FindDepoListFromTrack(self, trkId):
        x = [i for (i, depo) 
             in enumerate(self.depos) 
             if (trkId in depo.Contrib)]
            # if (trkId == depo.GetPrimaryId())]
        return x

    # ------------------------
    def PrintVertex(self):
        self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
        self.info['E_nu'] = self.genieTree.StdHepP4[3]*1000
        print(f"neutrino {self.info['nu_pdg']}: {self.info['E_nu']} MeV")
        posx = self.vertex.GetPosition().X()
        posy = self.vertex.GetPosition().Y()
        posz = self.vertex.GetPosition().Z()
        print(f"vertex @: ({posx}, {posy}, {posz}) [mm]")
        print(f"reaction: {self.vertex.GetReaction()}")
        # print(f"interaction #: {self.vertex.GetInteractionNumber()}")

        # print(f"{self.vertex.Particles.size()} particles at the vertex", )
        print(f'{"pdg":>8}{"name":>8}{"trkId":>6}{"mass":>10}{"KE":>10}')
        print(f'{"":>8}{"":>8}{"":>6}{"[MeV]":>10}{"[MeV]":>10}')
        print('-'*(8+8+6+10+10))
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            name = particle.GetName()
            trkId = particle.GetTrackId()
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass            
            print(f'{pdg:>8d}{name:>8s}{trkId:>6d}{mass:>10.2f}{KE:>10.2f}')
        print('-'*(8+8+6+10+10))

        print(f'{self.info}')

    # ------------------------
    def FillEnergyInfo(self):

        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            trkId = particle.GetTrackId()
            depoE = self.GetEnergyDepoWithDesendents(trkId)
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass            
            self.info['E_depoTotal'] += depoE
            # fill E_availList: lepton, proton, neutron, pi+-, pi0, others.
            if (pdg in [13, -13, 11, -11]):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][0] += (KE + mass)
                self.info['E_depoList'][0] += depoE
            elif (pdg == 2212):
                self.info['E_avail'] += KE
                self.info['E_availList'][1] += KE 
                self.info['E_depoList'][1] += depoE
            elif (pdg == 2112):
                self.info['E_avail'] += KE
                self.info['E_availList'][2] += KE
                self.info['E_depoList'][2] += depoE
            elif (pdg in [211, -211]):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][3] += (KE + mass)
                self.info['E_depoList'][3] += depoE
            elif (pdg == 111):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][4] += (KE + mass)
                self.info['E_depoList'][4] += depoE
            else:
                self.info['E_avail'] += KE
                self.info['E_availList'][5] += KE
                self.info['E_depoList'][5] += depoE

        # for track in self.tracks:
        #     pdg = track.GetPDGCode()
        #     depoE = track.energy['depoTotal']
        #     if (pdg in [13, -13, 11, -11]):
        #         self.info['E_depoList'][0] += depoE
        #     elif (pdg == 2212):
        #         self.info['E_depoList'][1] += depoE 
        #     elif (pdg == 2112):
        #         self.info['E_depoList'][2] += depoE
        #     elif (pdg in [211, -211]):
        #         self.info['E_depoList'][3] += depoE
        #     elif (pdg == 111):
        #         self.info['E_depoList'][4] += depoE
        #     else:
        #         self.info['E_depoList'][5] += depoE            


    # ------------------------
    def PrintDepo(self, i):
        # print(f"{self.depos.size} depo points stored in total")
        depo = self.depos[i]
        contrib = depo.GetContributors()
        primaryId = depo.GetPrimaryId()
        edep = depo.GetEnergyDeposit()
        # edepSecond = depo.GetSecondaryDeposit()
        trkLength = depo.GetTrackLength()
        print(contrib, "%6d %.2f %.2f" %(primaryId, edep, trkLength))

    # ------------------------
    def PrintTrack(self, trkId):
        self.PrintTracks(trkId, trkId+1)
        track = self.tracks[trkId]
        children = track.association['children']
        print(f'children: {children}')
        for childId in children:
            child = self.tracks[childId]
            name = child.GetName()
            mom = child.GetInitialMomentum()
            KE = mom.E() - mom.M()
            print(f"{childId} {name}: {KE:.2f} MeV, {child.energy['depoTotal']:.2f} MeV, {child.association['children']}")

        print(f"{track.Points.size()} points stored in track {trkId}")
        for point in track.Points:
            x = point.GetPosition().X()
            y = point.GetPosition().Y()
            z = point.GetPosition().Z() 
            t = point.GetPosition().T() 
            print(f"{point.GetProcess()}, {point.GetSubprocess()}, {x}, {y}, {z}, {t}")

    # ------------------------
    def PrintTracks(self, start=0, stop=-1):
        # print(f"{self.tracks.size} trajectories stored", )
        print(f"{'pdg':>8}{'name':>8}{'trkId':>6}{'parId':>6}{'acId':>6}{'KE':>10}{'selfDepo':>10}{'allDepo':>10}")
        print(f"{'':>8}{'':>8}{'':>6}{'':>6}{'':>6}{'[MeV]':>10}{'[MeV]':>10}{'[MeV]':>10}")
        print('-'*(8+8+6+6+6+10+10+10))

        for track in self.tracks[start:stop]:
            pdg = track.GetPDGCode()
            name = track.GetName()
            trkId = track.GetTrackId()
            parId = track.GetParentId()
            mom = track.GetInitialMomentum()
            mass = mom.M()
            KE = mom.E() - mass        
            ancestor = track.association['ancestor']
            selfDepo = track.energy['depoTotal']
            allDepo = self.GetEnergyDepoWithDesendents(trkId)
            print(f"{pdg:>8d}{name:>8s}{trkId:>6d}{parId:>6d}{ancestor:>6d}{KE:>10.2f}{selfDepo:>10.2f}{allDepo:>10.2f}")

        print('-'*(8+8+6+6+6+10+10+10))


    # ------------------------
    def Next(self):
        if self.currentEntry != self.nEntry -1:
            self.currentEntry += 1 
        else:
            self.currentEntry = 0
        self.Jump(self.currentEntry)

    # ------------------------
    def Prev(self):
        if self.currentEntry != 0:
            self.currentEntry -= 1 
        else:
            self.currentEntry = self.nEntry -1
        self.Jump(self.currentEntry)

    #-------------------------
    def GetEnergyDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        energy = track.energy['depoTotal']
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents(childId)
        return energy

    #-------------------------
    def GetEnergyDepoWithAncestor(self, acId):
        energy = 0
        for track in self.tracks:    
            ancestor = track.association['ancestor']
            if ancestor == acId:
                energy += track.energy['depoTotal']        
        return energy

    #-----------------------
    def GetFileName(self):
        return self.simTree.GetFile().GetName()
# ------------------------
if __name__ == "__main__":
    event = Event(sys.argv[1])
    event.Jump(0)
    event.PrintVertex()
    event.Next()
    event.PrintVertex()
    event.PrintTracks(0,8)    
    # event.PrintTrack(1)