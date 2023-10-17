import ROOT
from ROOT import TG4Event, TFile, TTree

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
        self.rootFile = TFile(self.fileName)
        self.simTree = self.rootFile.Get("EDepSimEvents")
        self.genieTree = self.rootFile.Get("DetSimPassThru/gRooTracker")

        self.nEntry = self.simTree.GetEntriesFast()
        if self.nEntry != self.genieTree.GetEntriesFast():
            print("Edep-sim tree and GENIE tree number of entries do not match!")
            sys.exit()        
        
        self.event = TG4Event()
        self.simTree.SetBranchAddress("Event", self.event)
    
    # ------------------------    
    def Jump(self, entryNo):
        print('reading event ', entryNo)

        self.currentEntry = entryNo
        self.simTree.GetEntry(entryNo)

        self.ReadVertex()
        self.ReadTracks()
        self.ReadEnergyDepo('SimEnergyDeposit')
    
    # ------------------------
    def ReadVertex(self):
        primaries = np.array(self.event.Primaries)
        if (primaries.size != 1):
            print("Number of primaries not equal to 1 (not neutrino vertex)!")
            return

        self.vertex = primaries[0]

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
            self.tracks[parId].association['children'].append(i)
            if parId == -1: continue
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

        E_tot = np.sum([depo.GetEnergyDeposit() for depo in self.depos])
        print('total deposit energy: ', E_tot)

    # ------------------------
    def FindDepoListFromTrack(self, trkId):
        x = [i for (i, depo) 
             in enumerate(self.depos) 
             if (trkId in depo.Contrib)]
            # if (trkId == depo.GetPrimaryId())]
        return x

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
        print(f"{track.Points.size()} points stored in track {trkId}")
        for point in track.Points:
            momx = point.GetMomentum().X()
            momy = point.GetMomentum().Y()
            momz = point.GetMomentum().Z() 
            print(f"{point.GetProcess()}, {point.GetSubprocess()}, {momx}, {momy}, {momz}")

    # ------------------------
    def PrintTracks(self, start=0, stop=-1):
        print(f"{self.tracks.size} trajectories stored", )
        print ("%8s%8s%6s%6s%10s%10s%10s" % ('pdg', 'name', 'trkId', 'parId', 'momx', 'momy', 'momz'))
        print ("%8s%8s%6s%6s%10s%10s%10s" % ('', '', '', '', '[MeV]', '[MeV]', '[MeV]'))
        print ("%8s%8s%6s%6s%10s%10s%10s" % ('-'*8, '-'*8, '-'*6, '-'*6, '-'*10, '-'*10, '-'*10))
        # if (stop == -1): stop = self.tracks.size
        for track in self.tracks[start:stop]:
            pdg = track.GetPDGCode()
            name = track.GetName()
            trkId = track.GetTrackId()
            parId = track.GetParentId()
            momx = track.GetInitialMomentum().X()
            momy = track.GetInitialMomentum().Y()
            momz = track.GetInitialMomentum().Z()
            # depoList = self.FindDepoListFromTrack(trkId)
            # depoEnergy = np.sum([depo.GetEnergyDeposit() for depo in self.depos[depoList]])    
            # if (depoEnergy<1): continue        
            print("%8d%8s%6d%6d%10.2f%10.2f%10.2f" % (pdg, name, trkId, parId, momx, momy, momz))
            # print(track.energy, track.association)
        print ("%8s%8s%6s%6s%10s%10s%10s" % ('-'*8, '-'*8, '-'*6, '-'*6, '-'*10, '-'*10, '-'*10))
    
    # ------------------------
    def PrintVertex(self):
        posx = self.vertex.GetPosition().X()
        posy = self.vertex.GetPosition().Y()
        posz = self.vertex.GetPosition().Z()
        print(f"vertex @: ({posx}, {posy}, {posz}) [mm]")
        print(f"reaction: {self.vertex.GetReaction()}")
        # print(f"interaction #: {self.vertex.GetInteractionNumber()}")
    
        # print(f"{self.vertex.Particles.size()} particles at the vertex", )
        print ("%8s%8s%6s%10s%10s" % ('pdg', 'name', 'trkId', 'mass', 'KE'))
        print ("%8s%8s%6s%10s%10s" % ('', '', '', '[MeV]', '[MeV]'))
        print ("%8s%8s%6s%10s%10s" % ('-'*8, '-'*8, '-'*6, '-'*10, '-'*10))
        for particle in self.vertex.Particles:
            pdg = particle.GetPDGCode()
            name = particle.GetName()
            trkId = particle.GetTrackId()
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass            
            print ("%8d%8s%6d%10.2f%10.2f" % (pdg, name, trkId, mass, KE))
        print ("%8s%8s%6s%10s%10s" % ('-'*8, '-'*8, '-'*6, '-'*10, '-'*10))

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

# ------------------------
if __name__ == "__main__":
    event = Event(sys.argv[1])
    event.Jump(0)
    event.PrintVertex()
    event.Next()
    event.PrintVertex()
    event.PrintTracks(0,8)    
    # event.PrintTrack(1)