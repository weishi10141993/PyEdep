import ROOT
from ROOT import TG4Event, TFile, TChain

import sys, os
import numpy as np
import random
import math

class Event:

    def __init__(self, fileName, evgen='Genie'):
        print("event: initilization")
        self.fileName = fileName
        self.evgen = evgen
        self.ReadTree()

        self.currentEntry = 0
        # create a folder to store plots
        self.plotpath = "./plots"
        if not os.path.exists( self.plotpath):
            os.makedirs( self.plotpath)
            print("plotpath '" + self.plotpath + "' did not exist. It has been created!")

    # ------------------------
    def ReadTree(self):
        # self.rootFile = TFile(self.fileName)
        self.simTree = TChain("EDepSimEvents")
        self.simTree.Add(self.fileName)
        self.nEntry = self.simTree.GetEntries()
        # Only Genie has gRooTracker, Marley doesn't
        if self.evgen == 'Genie':
            self.genieTree = TChain("DetSimPassThru/gRooTracker")
            self.genieTree.Add(self.fileName)
            # self.genieTree = self.rootFile.Get("DetSimPassThru/gRooTracker")
            if self.nEntry != self.genieTree.GetEntries():
                print("Edep-sim tree and GENIE tree number of entries do not match!")
                sys.exit()
        elif self.evgen == 'Marley':
            print("Marley events: Skip looking for gRooTracker directory.")
        else:
            print("Unknown event generator!")
            sys.exit()

        self.event = TG4Event()
        self.simTree.SetBranchAddress("Event", self.event)

    # ------------------------
    def Jump(self, entryNo):
        if entryNo%100 == 0:
            print(f'reading event {entryNo}/{self.nEntry}')

        self.currentEntry = entryNo
        self.simTree.GetEntry(entryNo)

        self.ReadVertex()
        self.ReadTracks()
        self.ReadEnergyDepo('SimEnergyDeposit')

        self.info = {}
        self.info['Event_ID'] = self.currentEntry
        self.info['E_nu'] = 0
        self.info['E_avail'] = 0
        self.info['E_availList'] = np.zeros(8) # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.info['E_depoTotal'] = 0
        self.info['E_depoTotal_track'] = 0
        self.info['Q_depoTotal_dots_th_75keV'] = 0
        self.info['E_depoTotal_no_ncap_depo'] = 0
        self.info['Q_depoTotal'] = 0
        self.info['Q_depoTotal_th_75keV'] = 0
        self.info['Q_depoTotal_th_500keV'] = 0
        self.info['Q_depoTotal_MBox'] = 0
        self.info['Q_depoTotal_MBox_th_75keV'] = 0
        self.info['Q_depoTotal_MBox_th_500keV'] = 0
        self.info['Q_depoTotal_no_ncap_depo'] = 0
        self.info['Q_depoTotal_th_75keV_no_ncap_depo'] = 0
        self.info['Q_depoTotal_th_500keV_no_ncap_depo'] = 0
        self.info['L_depoTotal_avg_220PEpMeV'] = 0
        self.info['L_depoTotal_avg_180PEpMeV'] = 0
        self.info['L_depoTotal_avg_140PEpMeV'] = 0
        self.info['L_depoTotal_avg_100PEpMeV'] = 0
        self.info['L_depoTotal_avg_35PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_220PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_180PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_140PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_100PEpMeV'] = 0
        self.info['L_depoTotal_MBox_avg_35PEpMeV'] = 0
        self.info['L_depoTotal_avg_220PEpMeV_no_ncap_depo'] = 0
        self.info['L_depoTotal_avg_180PEpMeV_no_ncap_depo'] = 0
        self.info['L_depoTotal_avg_140PEpMeV_no_ncap_depo'] = 0
        self.info['L_depoTotal_avg_100PEpMeV_no_ncap_depo'] = 0
        self.info['L_depoTotal_avg_35PEpMeV_no_ncap_depo'] = 0
        self.info['E_depoList'] = np.zeros(8) # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.info['E_depoList_track'] = np.zeros(8)
        self.info['E_depoList_no_ncap_depo'] = np.zeros(8)
        self.info['Q_depoList']               = np.zeros(8)
        self.info['Q_depoList_th_75keV']      = np.zeros(8)
        self.info['Q_depoList_th_500keV']     = np.zeros(8)
        self.info['Q_depoList_dots_th_75keV'] = np.zeros(8)
        self.info['Q_depoList_MBox']           = np.zeros(8)
        self.info['Q_depoList_MBox_th_75keV']  = np.zeros(8)
        self.info['Q_depoList_MBox_th_500keV'] = np.zeros(8)
        self.info['Q_depoList_no_ncap_depo']               = np.zeros(8)
        self.info['Q_depoList_th_75keV_no_ncap_depo']      = np.zeros(8)
        self.info['Q_depoList_th_500keV_no_ncap_depo']     = np.zeros(8)
        self.info['L_depoList_avg_220PEpMeV'] = np.zeros(8)
        self.info['L_depoList_avg_180PEpMeV'] = np.zeros(8)
        self.info['L_depoList_avg_140PEpMeV'] = np.zeros(8)
        self.info['L_depoList_avg_100PEpMeV'] = np.zeros(8)
        self.info['L_depoList_avg_35PEpMeV']  = np.zeros(8)
        self.info['L_depoList_MBox_avg_220PEpMeV'] = np.zeros(8)
        self.info['L_depoList_MBox_avg_180PEpMeV'] = np.zeros(8)
        self.info['L_depoList_MBox_avg_140PEpMeV'] = np.zeros(8)
        self.info['L_depoList_MBox_avg_100PEpMeV'] = np.zeros(8)
        self.info['L_depoList_MBox_avg_35PEpMeV']  = np.zeros(8)
        self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'] = np.zeros(8)
        self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'] = np.zeros(8)
        self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'] = np.zeros(8)
        self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'] = np.zeros(8)
        self.info['L_depoList_avg_35PEpMeV_no_ncap_depo']  = np.zeros(8)
        self.info['N_parList'] = np.zeros(8) # lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
        self.info['nu_pdg'] = 0
        self.info['nu_xs'] = self.vertex.GetCrossSection()
        self.info['nu_proc'], self.info['nu_nucl'] = self.GetReaction()
        self.FillEnergyInfo()
        if self.evgen == 'Genie':
            self.ReadGenie()
        elif self.evgen == 'Marley':
            self.ReadMarley()
        else:
            print("Unknown event generator!")
            sys.exit()

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
    def ReadMarley(self):
        # we are mostly looking at nue anyway
        self.info['nu_pdg'] = 12
        if self.GetnuPDGFromFileName() == 'nue': self.info['nu_pdg'] = 12
        if self.GetnuPDGFromFileName() == 'numu': self.info['nu_pdg'] = 14
        if self.GetnuPDGFromFileName() == 'anue': self.info['nu_pdg'] = -12
        if self.GetnuPDGFromFileName() == 'anumu': self.info['nu_pdg'] = -14
        self.info['E_nu'] = self.GetEnuFromFileName()
        if self.currentEntry == 0: # only print once
            print("Marley events: assert info from file name")

    # ------------------------
    def ReadVertex(self):
        primaries = np.array(self.event.Primaries)
        if (primaries.size != 1 and self.evgen == 'Genie'):
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
            self.tracks[i].length = {}
            self.tracks[i].energy['depoTotal'] = 0
            self.tracks[i].energy['depoTotal_no_ncap_depo'] = 0
            self.tracks[i].energy['depoTotal_charge'] = 0
            self.tracks[i].energy['depoTotal_charge_MBox'] = 0
            self.tracks[i].energy['depoTotal_charge_no_ncap_depo'] = 0
            self.tracks[i].energy['depoTotal_charge_th_75keV'] = 0
            self.tracks[i].energy['depoTotal_charge_MBox_th_75keV'] = 0
            self.tracks[i].energy['depoTotal_charge_th_75keV_no_ncap_depo'] = 0
            self.tracks[i].energy['depoTotal_charge_th_500keV'] = 0
            self.tracks[i].energy['depoTotal_charge_MBox_th_500keV'] = 0
            self.tracks[i].energy['depoTotal_charge_th_500keV_no_ncap_depo'] = 0
            self.tracks[i].energy['depoTotal_light_avg_220PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_220PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_220PEpMeV_exclude_ncap_edepo'] = 0
            self.tracks[i].energy['depoTotal_light_avg_180PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_180PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_180PEpMeV_exclude_ncap_edepo'] = 0
            self.tracks[i].energy['depoTotal_light_avg_140PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_140PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_140PEpMeV_exclude_ncap_edepo'] = 0
            self.tracks[i].energy['depoTotal_light_avg_100PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_100PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_100PEpMeV_exclude_ncap_edepo'] = 0
            self.tracks[i].energy['depoTotal_light_avg_35PEpMeV'] = 0
            self.tracks[i].energy['depoTotal_light_avg_35PEpMeV_MBox'] = 0
            self.tracks[i].energy['depoTotal_light_avg_35PEpMeV_exclude_ncap_edepo'] = 0
            self.tracks[i].length['selfDepo'] = 0
            #self.tracks[i].length['selfDepo_no_ncap_depo'] = 0
            self.tracks[i].association['depoList'] = []
            #self.tracks[i].association['depoList_no_ncap_depo'] = []
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
        mm2cm = 0.1

        for i, depo in enumerate(self.depos):
            trkId = depo.Contrib[0]
            edep = depo.GetEnergyDeposit()
            trkLength = depo.GetTrackLength() *mm2cm
            Qdep = self.ChargeBirksLaw(edep, trkLength)
            Qdep_MBox = self.ChargeModifiedBoxModel(edep, trkLength) # Modified Box Model
            track = self.tracks[trkId]
            track.association['depoList'].append(i)
            track.energy['depoTotal'] += edep
            track.energy['depoTotal_charge'] += Qdep
            track.energy['depoTotal_charge_MBox'] += Qdep_MBox
            track.length['selfDepo'] += trkLength
            if Qdep > 0.075: track.energy['depoTotal_charge_th_75keV'] += Qdep
            if Qdep > 0.5: track.energy['depoTotal_charge_th_500keV'] += Qdep
            if Qdep_MBox > 0.075: track.energy['depoTotal_charge_MBox_th_75keV'] += Qdep_MBox
            if Qdep_MBox > 0.5: track.energy['depoTotal_charge_MBox_th_500keV'] += Qdep_MBox

            # Light deposit
            # edepSecond = depo.GetSecondaryDeposit()
            # Ideally this should be the light dep
            # how the light yield is actually derived: 25k photons / MeV are simualted and 220 on average are detected
            # the overall photon collection efficiency (PCE) is
            #    220/25000=0.88% PCE
            #    180/25000=0.72% PCE
            #    140/25000=0.56% PCE
            #    100/25000=0.4% PCE
            #    35/25000=0.14% PCE

            #    220/21622=1.0% PCE
            #    180/21622=0.83% PCE
            #    140/21622=0.65% PCE
            #    100/21622=0.46% PCE
            #    35/21622=0.16% PCE
            # (dL/W_ph)*PCE = (dL/19.5eV)*PCE = number of photons (taking into account Birks model)
            Ldep = edep - Qdep
            Ldep_MBox = edep - Qdep_MBox

            nPE_220 = (Ldep*1000000/19.5)*0.01
            nPE_180 = (Ldep*1000000/19.5)*0.0083
            nPE_140 = (Ldep*1000000/19.5)*0.0065
            nPE_100 = (Ldep*1000000/19.5)*0.0046
            nPE_35  = (Ldep*1000000/19.5)*0.0016

            nPE_220_MBox = Ldep_MBox*220
            nPE_180_MBox = Ldep_MBox*180
            nPE_140_MBox = Ldep_MBox*140
            nPE_100_MBox = Ldep_MBox*100
            nPE_35_MBox  = Ldep_MBox*35

            nPE_220_detected = random.gauss(nPE_220, math.sqrt(nPE_220))
            nPE_180_detected = random.gauss(nPE_180, math.sqrt(nPE_180))
            nPE_140_detected = random.gauss(nPE_140, math.sqrt(nPE_140))
            nPE_100_detected = random.gauss(nPE_100, math.sqrt(nPE_100))
            nPE_35_detected  = random.gauss(nPE_35, math.sqrt(nPE_35))

            nPE_220_MBox_detected = random.gauss(nPE_220_MBox, math.sqrt(nPE_220_MBox))
            nPE_180_MBox_detected = random.gauss(nPE_180_MBox, math.sqrt(nPE_180_MBox))
            nPE_140_MBox_detected = random.gauss(nPE_140_MBox, math.sqrt(nPE_140_MBox))
            nPE_100_MBox_detected = random.gauss(nPE_100_MBox, math.sqrt(nPE_100_MBox))
            nPE_35_MBox_detected  = random.gauss(nPE_35_MBox,  math.sqrt(nPE_35_MBox))

            Ldep_avg_220PEpMeV_detected = nPE_220_detected*19.5/1000000/0.01 # MeV
            Ldep_avg_180PEpMeV_detected = nPE_180_detected*19.5/1000000/0.0083
            Ldep_avg_140PEpMeV_detected = nPE_140_detected*19.5/1000000/0.0065
            Ldep_avg_100PEpMeV_detected = nPE_100_detected*19.5/1000000/0.0046
            Ldep_avg_35PEpMeV_detected  = nPE_35_detected*19.5/1000000/0.0016

            Ldep_avg_220PEpMeV_MBox_detected = nPE_220_MBox_detected/220
            Ldep_avg_180PEpMeV_MBox_detected = nPE_180_MBox_detected/180
            Ldep_avg_140PEpMeV_MBox_detected = nPE_140_MBox_detected/140
            Ldep_avg_100PEpMeV_MBox_detected = nPE_100_MBox_detected/100
            Ldep_avg_35PEpMeV_MBox_detected  = nPE_35_MBox_detected/35

            track.energy['depoTotal_light_avg_220PEpMeV'] += Ldep_avg_220PEpMeV_detected
            track.energy['depoTotal_light_avg_180PEpMeV'] += Ldep_avg_180PEpMeV_detected
            track.energy['depoTotal_light_avg_140PEpMeV'] += Ldep_avg_140PEpMeV_detected
            track.energy['depoTotal_light_avg_100PEpMeV'] += Ldep_avg_100PEpMeV_detected
            track.energy['depoTotal_light_avg_35PEpMeV']  += Ldep_avg_35PEpMeV_detected

            track.energy['depoTotal_light_avg_220PEpMeV_MBox'] += Ldep_avg_220PEpMeV_MBox_detected
            track.energy['depoTotal_light_avg_180PEpMeV_MBox'] += Ldep_avg_180PEpMeV_MBox_detected
            track.energy['depoTotal_light_avg_140PEpMeV_MBox'] += Ldep_avg_140PEpMeV_MBox_detected
            track.energy['depoTotal_light_avg_100PEpMeV_MBox'] += Ldep_avg_100PEpMeV_MBox_detected
            track.energy['depoTotal_light_avg_35PEpMeV_MBox']  += Ldep_avg_35PEpMeV_MBox_detected

            ####################################
            # Exclude neutron capture e deposits
            ####################################
            if self.EdepoFromNeutronCap(trkId):
                # skip this edepo for Q and L calculation
                continue
                # if the edepo is from neutron cap, then set this track e to zero, can't detect
                #track.energy['depoTotal_light_avg_220PEpMeV'] = 0
                #track.energy['depoTotal_light_avg_180PEpMeV'] = 0
                #track.energy['depoTotal_light_avg_140PEpMeV'] = 0
                #track.energy['depoTotal_light_avg_100PEpMeV'] = 0
                #track.energy['depoTotal_light_avg_35PEpMeV']  = 0

                #track.energy['depoTotal'] = 0
                #track.energy['depoTotal_charge'] = 0
                #if Qdep > 0.075: track.energy['depoTotal_charge_th_75keV'] = 0
                #if Qdep > 0.5: track.energy['depoTotal_charge_th_500keV'] = 0


            track.energy['depoTotal_light_avg_220PEpMeV_exclude_ncap_edepo'] += Ldep_avg_220PEpMeV_detected
            track.energy['depoTotal_light_avg_180PEpMeV_exclude_ncap_edepo'] += Ldep_avg_180PEpMeV_detected
            track.energy['depoTotal_light_avg_140PEpMeV_exclude_ncap_edepo'] += Ldep_avg_140PEpMeV_detected
            track.energy['depoTotal_light_avg_100PEpMeV_exclude_ncap_edepo'] += Ldep_avg_100PEpMeV_detected
            track.energy['depoTotal_light_avg_35PEpMeV_exclude_ncap_edepo']  += Ldep_avg_35PEpMeV_detected

            track.energy['depoTotal_no_ncap_depo'] += edep
            track.energy['depoTotal_charge_no_ncap_depo'] += Qdep
            if Qdep > 0.075: track.energy['depoTotal_charge_th_75keV_no_ncap_depo'] += Qdep
            if Qdep > 0.5: track.energy['depoTotal_charge_th_500keV_no_ncap_depo'] += Qdep

        # E_tot = np.sum([depo.GetEnergyDeposit() for depo in self.depos])
        # print('total deposit energy: ', E_tot)

    # ------------------------
    def ChargeBirksLaw(self, edep, trkLength):
        # PHYS. REV. D 99, 036009 (2019)
        # return R = (dQ/dx)/(dE/dx) = dQ/dE = A/(1+kQ*dE/dx)
        # A = 0.8, kQ = 0.0972 g/(MeV*cm2) @500V/cm
        # LAr density 1.4g/cm^3
        # when considering light, need to take account excitation and ionization ratio (0.83 vs 0.17)
        return 0.83*edep*0.8/(1+0.0972*edep/trkLength/1.4)

    # ------------------------
    def ChargeModifiedBoxModel(self, edep, trkLength):
        # https://lar.bnl.gov/properties/pass.html
        # return R = (dQ/dx)/(dE/dx) = dQ/dE = ln(A+B*dE/dx)/(B*dE/dx)
        # A = 0.930, B = 0.424 g/(MeV*cm2) @500V/cm
        # LAr density 1.4g/cm^3
        # when considering light, need to take account excitation and ionization ratio (0.83 vs 0.17)
        return 0.83*edep*np.log(0.930+0.424*edep/trkLength/1.4)/(0.424*edep/trkLength/1.4)

    # ------------------------
    def FindDepoListFromTrack(self, trkId):
        x = [i for (i, depo)
             in enumerate(self.depos)
             if (trkId in depo.Contrib)]
            # if (trkId == depo.GetPrimaryId())]
        return x

    # ------------------------
    def EdepoFromNeutronCap(self, trkId):
        # if edepo associated track's parent trk belong to/can trace back to neutron
        # and the neutron has a last process and subprocess 4-131 (i.e. G4 n capture), skip the edepo
        track = self.tracks[trkId]
        parId = track.GetParentId()
        partrack = self.tracks[parId]
        parpdg = partrack.GetPDGCode()
        flag = False

        if parpdg == 2112:
            # check this neutron's last process and subprocess
            if partrack.Points[partrack.Points.size()-1].GetProcess() == 4 and partrack.Points[partrack.Points.size()-1].GetSubprocess() == 131:
                #print(f"{partrack.Points[partrack.Points.size()-1].GetProcess()}, {partrack.Points[partrack.Points.size()-1].GetSubprocess()}")
                flag = True
            else:
                flag = False
        elif parId == -1:
            # hit the beginning
            flag = False
        else:
            # if not neutron keep tracing backward
            flag = self.EdepoFromNeutronCap(parId)

        return flag

    # ------------------------
    def PrintVertex(self):
        if self.evgen == 'Genie':
            self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
            self.info['E_nu'] = self.genieTree.StdHepP4[3]*1000
            print(f"neutrino {self.info['nu_pdg']}: {self.info['E_nu']} MeV")
        elif self.evgen == 'Marley':
            self.ReadMarley()
        else:
            print("Unknown event generator!")
            sys.exit()

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
            trkId = particle.GetTrackId()
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            if trkId < 0: continue
            pdg = particle.GetPDGCode()
            name = particle.GetName()
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            print(f'{pdg:>8d}{name:>8s}{trkId:>6d}{mass:>10.2f}{KE:>10.2f}')
        print('-'*(8+8+6+10+10))

        print(f'{self.info}')

    # ------------------------
    def FillEnergyInfo(self):

        for particle in self.vertex.Particles:
            trkId = particle.GetTrackId()
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            if trkId < 0: continue
            pdg = particle.GetPDGCode()
            depoE = self.GetEnergyDepoWithDesendents(trkId)
            depoE_no_ncap_depo = self.GetEnergyDepoWithDesendentsExcludeNCapEDepo(trkId)
            track = self.tracks[trkId]
            # if longer than 2cm, assume can reconstruct dE
            if track.length['selfDepo'] > 2:
                depoE_track = track.energy['depoTotal']
                depoQ_dots = 0
            else:
                # dots/blips
                depoE_track = 0
                depoQ_dots = track.energy['depoTotal_charge_th_75keV']
            depoQ        = self.GetChargeDepoWithDesendents(trkId)[0]
            depoQ_75keV  = self.GetChargeDepoWithDesendents(trkId)[1]
            depoQ_500keV = self.GetChargeDepoWithDesendents(trkId)[2]
            depoQ_MBox        = self.GetChargeMBoxDepoWithDesendents(trkId)[0]
            depoQ_MBox_75keV  = self.GetChargeMBoxDepoWithDesendents(trkId)[1]
            depoQ_MBox_500keV = self.GetChargeMBoxDepoWithDesendents(trkId)[2]
            depoQ_no_ncap_depo        = self.GetChargeDepoWithDesendentsExcludeNCapEDepo(trkId)[0]
            depoQ_75keV_no_ncap_depo  = self.GetChargeDepoWithDesendentsExcludeNCapEDepo(trkId)[1]
            depoQ_500keV_no_ncap_depo = self.GetChargeDepoWithDesendentsExcludeNCapEDepo(trkId)[2]

            depoL_220PEpMeV = self.GetLightDepoWithDesendents(trkId)[0]
            depoL_180PEpMeV = self.GetLightDepoWithDesendents(trkId)[1]
            depoL_140PEpMeV = self.GetLightDepoWithDesendents(trkId)[2]
            depoL_100PEpMeV = self.GetLightDepoWithDesendents(trkId)[3]
            depoL_35PEpMeV  = self.GetLightDepoWithDesendents(trkId)[4]
            depoL_MBox_220PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[0]
            depoL_MBox_180PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[1]
            depoL_MBox_140PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[2]
            depoL_MBox_100PEpMeV = self.GetLightMBoxDepoWithDesendents(trkId)[3]
            depoL_MBox_35PEpMeV  = self.GetLightMBoxDepoWithDesendents(trkId)[4]
            depoL_220PEpMeV_exclude_ncap_depo = self.GetLightDepoWithDesendentsExcludeNCapEDepo(trkId)[0]
            depoL_180PEpMeV_exclude_ncap_depo = self.GetLightDepoWithDesendentsExcludeNCapEDepo(trkId)[1]
            depoL_140PEpMeV_exclude_ncap_depo = self.GetLightDepoWithDesendentsExcludeNCapEDepo(trkId)[2]
            depoL_100PEpMeV_exclude_ncap_depo = self.GetLightDepoWithDesendentsExcludeNCapEDepo(trkId)[3]
            depoL_35PEpMeV_exclude_ncap_depo  = self.GetLightDepoWithDesendentsExcludeNCapEDepo(trkId)[4]

            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            self.info['E_depoTotal'] += depoE
            self.info['E_depoTotal_track'] += depoE_track
            self.info['E_depoTotal_no_ncap_depo'] += depoE_no_ncap_depo
            self.info['Q_depoTotal'] += depoQ
            self.info['Q_depoTotal_th_75keV'] += depoQ_75keV
            self.info['Q_depoTotal_th_500keV'] += depoQ_500keV
            self.info['Q_depoTotal_dots_th_75keV'] += depoQ_dots
            self.info['Q_depoTotal_MBox'] += depoQ_MBox
            self.info['Q_depoTotal_MBox_th_75keV'] += depoQ_MBox_75keV
            self.info['Q_depoTotal_MBox_th_500keV'] += depoQ_MBox_500keV
            self.info['Q_depoTotal_no_ncap_depo'] += depoQ_no_ncap_depo
            self.info['Q_depoTotal_th_75keV_no_ncap_depo'] += depoQ_75keV_no_ncap_depo
            self.info['Q_depoTotal_th_500keV_no_ncap_depo'] += depoQ_500keV_no_ncap_depo

            self.info['L_depoTotal_avg_220PEpMeV'] += depoL_220PEpMeV
            self.info['L_depoTotal_avg_180PEpMeV'] += depoL_180PEpMeV
            self.info['L_depoTotal_avg_140PEpMeV'] += depoL_140PEpMeV
            self.info['L_depoTotal_avg_100PEpMeV'] += depoL_100PEpMeV
            self.info['L_depoTotal_avg_35PEpMeV']  += depoL_35PEpMeV
            self.info['L_depoTotal_MBox_avg_220PEpMeV'] += depoL_MBox_220PEpMeV
            self.info['L_depoTotal_MBox_avg_180PEpMeV'] += depoL_MBox_180PEpMeV
            self.info['L_depoTotal_MBox_avg_140PEpMeV'] += depoL_MBox_140PEpMeV
            self.info['L_depoTotal_MBox_avg_100PEpMeV'] += depoL_MBox_100PEpMeV
            self.info['L_depoTotal_MBox_avg_35PEpMeV']  += depoL_MBox_35PEpMeV
            self.info['L_depoTotal_avg_220PEpMeV_no_ncap_depo'] += depoL_220PEpMeV_exclude_ncap_depo
            self.info['L_depoTotal_avg_180PEpMeV_no_ncap_depo'] += depoL_180PEpMeV_exclude_ncap_depo
            self.info['L_depoTotal_avg_140PEpMeV_no_ncap_depo'] += depoL_140PEpMeV_exclude_ncap_depo
            self.info['L_depoTotal_avg_100PEpMeV_no_ncap_depo'] += depoL_100PEpMeV_exclude_ncap_depo
            self.info['L_depoTotal_avg_35PEpMeV_no_ncap_depo']  += depoL_35PEpMeV_exclude_ncap_depo
            # fill E_availList: lepton, proton, neutron, pi+-, pi0, gamma, alpha, others.
            if (pdg in [13, -13, 11, -11]):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][0] += (KE + mass)
                self.info['E_depoList'][0]       += depoE
                self.info['E_depoList_track'][0] += depoE_track
                self.info['E_depoList_no_ncap_depo'][0]   += depoE_no_ncap_depo
                self.info['Q_depoList'][0]                += depoQ
                self.info['Q_depoList_th_75keV'][0]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][0]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][0]  += depoQ_dots
                self.info['Q_depoList_MBox'][0]           += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][0]  += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][0] += depoQ_MBox_500keV
                self.info['Q_depoList_no_ncap_depo'][0]                += depoQ_no_ncap_depo
                self.info['Q_depoList_th_75keV_no_ncap_depo'][0]       += depoQ_75keV_no_ncap_depo
                self.info['Q_depoList_th_500keV_no_ncap_depo'][0]      += depoQ_500keV_no_ncap_depo
                self.info['L_depoList_avg_220PEpMeV'][0] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][0] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][0] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][0] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][0]  += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][0] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][0] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][0] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][0] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][0]  += depoL_MBox_35PEpMeV
                self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'][0] += depoL_220PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'][0] += depoL_180PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'][0] += depoL_140PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'][0] += depoL_100PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_35PEpMeV_no_ncap_depo'][0]  += depoL_35PEpMeV_exclude_ncap_depo
                self.info['N_parList'][0] += 1
            elif (pdg == 2212): # proton
                self.info['E_avail'] += KE
                self.info['E_availList'][1] += KE
                self.info['E_depoList'][1]       += depoE
                self.info['E_depoList_track'][1] += depoE_track
                self.info['E_depoList_no_ncap_depo'][1]   += depoE_no_ncap_depo
                self.info['Q_depoList'][1]                += depoQ
                self.info['Q_depoList_th_75keV'][1]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][1]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][1]  += depoQ_dots
                self.info['Q_depoList_MBox'][1]           += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][1]  += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][1] += depoQ_MBox_500keV
                self.info['Q_depoList_no_ncap_depo'][1]                += depoQ_no_ncap_depo
                self.info['Q_depoList_th_75keV_no_ncap_depo'][1]       += depoQ_75keV_no_ncap_depo
                self.info['Q_depoList_th_500keV_no_ncap_depo'][1]      += depoQ_500keV_no_ncap_depo
                self.info['L_depoList_avg_220PEpMeV'][1] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][1] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][1] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][1] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][1]  += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][1] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][1] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][1] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][1] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][1]  += depoL_MBox_35PEpMeV
                self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'][1] += depoL_220PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'][1] += depoL_180PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'][1] += depoL_140PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'][1] += depoL_100PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_35PEpMeV_no_ncap_depo'][1]  += depoL_35PEpMeV_exclude_ncap_depo
                self.info['N_parList'][1] += 1
            elif (pdg == 2112): # neutron
                self.info['E_avail'] += KE
                self.info['E_availList'][2] += KE
                self.info['E_depoList'][2]       += depoE
                self.info['E_depoList_track'][2] += depoE_track
                self.info['E_depoList_no_ncap_depo'][2]   += depoE_no_ncap_depo
                self.info['Q_depoList'][2]                += depoQ
                self.info['Q_depoList_th_75keV'][2]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][2]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][2]  += depoQ_dots
                self.info['Q_depoList_MBox'][2]           += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][2]  += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][2] += depoQ_MBox_500keV
                self.info['Q_depoList_no_ncap_depo'][2]                += depoQ_no_ncap_depo
                self.info['Q_depoList_th_75keV_no_ncap_depo'][2]       += depoQ_75keV_no_ncap_depo
                self.info['Q_depoList_th_500keV_no_ncap_depo'][2]      += depoQ_500keV_no_ncap_depo
                self.info['L_depoList_avg_220PEpMeV'][2] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][2] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][2] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][2] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][2]  += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][2] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][2] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][2] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][2] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][2]  += depoL_MBox_35PEpMeV
                self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'][2] += depoL_220PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'][2] += depoL_180PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'][2] += depoL_140PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'][2] += depoL_100PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_35PEpMeV_no_ncap_depo'][2]  += depoL_35PEpMeV_exclude_ncap_depo
                self.info['N_parList'][2] += 1
            elif (pdg in [211, -211]): #pi+/-
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][3] += (KE + mass)
                self.info['E_depoList'][3]       += depoE
                self.info['E_depoList_track'][3] += depoE_track
                self.info['E_depoList_no_ncap_depo'][3]   += depoE_no_ncap_depo
                self.info['Q_depoList'][3]                += depoQ
                self.info['Q_depoList_th_75keV'][3]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][3]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][3]  += depoQ_dots
                self.info['Q_depoList_MBox'][3]           += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][3]  += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][3] += depoQ_MBox_500keV
                self.info['Q_depoList_no_ncap_depo'][3]                += depoQ_no_ncap_depo
                self.info['Q_depoList_th_75keV_no_ncap_depo'][3]       += depoQ_75keV_no_ncap_depo
                self.info['Q_depoList_th_500keV_no_ncap_depo'][3]      += depoQ_500keV_no_ncap_depo
                self.info['L_depoList_avg_220PEpMeV'][3] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][3] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][3] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][3] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][3]  += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][3] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][3] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][3] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][3] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][3]  += depoL_MBox_35PEpMeV
                self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'][3] += depoL_220PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'][3] += depoL_180PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'][3] += depoL_140PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'][3] += depoL_100PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_35PEpMeV_no_ncap_depo'][3]  += depoL_35PEpMeV_exclude_ncap_depo
                self.info['N_parList'][3] += 1
            elif (pdg == 111): #pi0
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][4] += (KE + mass)
                self.info['E_depoList'][4]       += depoE
                self.info['E_depoList_track'][4] += depoE_track
                self.info['E_depoList_no_ncap_depo'][4]   += depoE_no_ncap_depo
                self.info['Q_depoList'][4]                += depoQ
                self.info['Q_depoList_th_75keV'][4]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][4]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][4]  += depoQ_dots
                self.info['Q_depoList_MBox'][4]           += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][4]  += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][4] += depoQ_MBox_500keV
                self.info['Q_depoList_no_ncap_depo'][4]                += depoQ_no_ncap_depo
                self.info['Q_depoList_th_75keV_no_ncap_depo'][4]       += depoQ_75keV_no_ncap_depo
                self.info['Q_depoList_th_500keV_no_ncap_depo'][4]      += depoQ_500keV_no_ncap_depo
                self.info['L_depoList_avg_220PEpMeV'][4] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][4] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][4] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][4] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][4]  += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][4] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][4] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][4] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][4] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][4]  += depoL_MBox_35PEpMeV
                self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'][4] += depoL_220PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'][4] += depoL_180PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'][4] += depoL_140PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'][4] += depoL_100PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_35PEpMeV_no_ncap_depo'][4]  += depoL_35PEpMeV_exclude_ncap_depo
                self.info['N_parList'][4] += 1
            elif (pdg == 22): #gamma
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][5] += (KE + mass)
                self.info['E_depoList'][5]       += depoE
                self.info['E_depoList_track'][5] += depoE_track
                self.info['E_depoList_no_ncap_depo'][5]   += depoE_no_ncap_depo
                self.info['Q_depoList'][5]                += depoQ
                self.info['Q_depoList_th_75keV'][5]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][5]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][5]  += depoQ_dots
                self.info['Q_depoList_MBox'][5]           += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][5]  += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][5] += depoQ_MBox_500keV
                self.info['Q_depoList_no_ncap_depo'][5]                += depoQ_no_ncap_depo
                self.info['Q_depoList_th_75keV_no_ncap_depo'][5]       += depoQ_75keV_no_ncap_depo
                self.info['Q_depoList_th_500keV_no_ncap_depo'][5]      += depoQ_500keV_no_ncap_depo
                self.info['L_depoList_avg_220PEpMeV'][5] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][5] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][5] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][5] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][5] += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][5] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][5] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][5] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][5] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][5]  += depoL_MBox_35PEpMeV
                self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'][5] += depoL_220PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'][5] += depoL_180PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'][5] += depoL_140PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'][5] += depoL_100PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_35PEpMeV_no_ncap_depo'][5]  += depoL_35PEpMeV_exclude_ncap_depo
                self.info['N_parList'][5] += 1
            elif (pdg == 1000020040): #alpha
                self.info['E_avail'] += KE
                self.info['E_availList'][6] += KE
                self.info['E_depoList'][6]       += depoE
                self.info['E_depoList_track'][6] += depoE_track
                self.info['E_depoList_no_ncap_depo'][6]   += depoE_no_ncap_depo
                self.info['Q_depoList'][6]                += depoQ
                self.info['Q_depoList_th_75keV'][6]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][6]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][6]  += depoQ_dots
                self.info['Q_depoList_MBox'][6]           += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][6]  += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][6] += depoQ_MBox_500keV
                self.info['Q_depoList_no_ncap_depo'][6]                += depoQ_no_ncap_depo
                self.info['Q_depoList_th_75keV_no_ncap_depo'][6]       += depoQ_75keV_no_ncap_depo
                self.info['Q_depoList_th_500keV_no_ncap_depo'][6]      += depoQ_500keV_no_ncap_depo
                self.info['L_depoList_avg_220PEpMeV'][6] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][6] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][6] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][6] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][6]  += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][6] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][6] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][6] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][6] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][6]  += depoL_MBox_35PEpMeV
                self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'][6] += depoL_220PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'][6] += depoL_180PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'][6] += depoL_140PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'][6] += depoL_100PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_35PEpMeV_no_ncap_depo'][6]  += depoL_35PEpMeV_exclude_ncap_depo
                self.info['N_parList'][6] += 1
            else:
                self.info['E_avail'] += KE
                self.info['E_availList'][7] += KE
                self.info['E_depoList'][7]       += depoE
                self.info['E_depoList_track'][7] += depoE_track
                self.info['E_depoList_no_ncap_depo'][7]   += depoE_no_ncap_depo
                self.info['Q_depoList'][7]                += depoQ
                self.info['Q_depoList_th_75keV'][7]       += depoQ_75keV
                self.info['Q_depoList_th_500keV'][7]      += depoQ_500keV
                self.info['Q_depoList_dots_th_75keV'][7]  += depoQ_dots
                self.info['Q_depoList_MBox'][7]           += depoQ_MBox
                self.info['Q_depoList_MBox_th_75keV'][7]  += depoQ_MBox_75keV
                self.info['Q_depoList_MBox_th_500keV'][7] += depoQ_MBox_500keV
                self.info['Q_depoList_no_ncap_depo'][7]                += depoQ_no_ncap_depo
                self.info['Q_depoList_th_75keV_no_ncap_depo'][7]       += depoQ_75keV_no_ncap_depo
                self.info['Q_depoList_th_500keV_no_ncap_depo'][7]      += depoQ_500keV_no_ncap_depo
                self.info['L_depoList_avg_220PEpMeV'][7] += depoL_220PEpMeV
                self.info['L_depoList_avg_180PEpMeV'][7] += depoL_180PEpMeV
                self.info['L_depoList_avg_140PEpMeV'][7] += depoL_140PEpMeV
                self.info['L_depoList_avg_100PEpMeV'][7] += depoL_100PEpMeV
                self.info['L_depoList_avg_35PEpMeV'][7]  += depoL_35PEpMeV
                self.info['L_depoList_MBox_avg_220PEpMeV'][7] += depoL_MBox_220PEpMeV
                self.info['L_depoList_MBox_avg_180PEpMeV'][7] += depoL_MBox_180PEpMeV
                self.info['L_depoList_MBox_avg_140PEpMeV'][7] += depoL_MBox_140PEpMeV
                self.info['L_depoList_MBox_avg_100PEpMeV'][7] += depoL_MBox_100PEpMeV
                self.info['L_depoList_MBox_avg_35PEpMeV'][7]  += depoL_MBox_35PEpMeV
                self.info['L_depoList_avg_220PEpMeV_no_ncap_depo'][7] += depoL_220PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_180PEpMeV_no_ncap_depo'][7] += depoL_180PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_140PEpMeV_no_ncap_depo'][7] += depoL_140PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_100PEpMeV_no_ncap_depo'][7] += depoL_100PEpMeV_exclude_ncap_depo
                self.info['L_depoList_avg_35PEpMeV_no_ncap_depo'][7]  += depoL_35PEpMeV_exclude_ncap_depo
                self.info['N_parList'][7] += 1

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
        #print(f'children: {children}')
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

        #print(f"{point.GetProcess()}, {point.GetSubprocess()}, {x}, {y}, {z}, {t}")
        return t # last point is capture time

    # ------------------------
    def PrintTracks(self, start=0, stop=-1):
        # print(f"{self.tracks.size} trajectories stored", )

        print(f"{'pdg':>8}{'name':>8}{'trkId':>6}{'parId':>6}{'acId':>6}{'KE':>10}{'selfDepo':>10}{'allDepo':>10}")
        print(f"{'':>8}{'':>8}{'':>6}{'':>6}{'':>6}{'[MeV]':>10}{'[MeV]':>10}{'[MeV]':>10}")
        print('-'*(8+8+6+6+6+10+10+10))

        neutrontrkId = -1
        neutronKE = -1

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
            # Take the last neutron and check its KE and capture time
            if pdg == 2112:
                neutrontrkId = trkId
                neutronKE = KE

        #print('-'*(8+8+6+6+6+10+10+10))
        return [neutrontrkId, neutronKE]


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
    def GetEnergyDepoWithDesendentsExcludeNCapEDepo(self, trkId):
        track = self.tracks[trkId]
        energy_no_ncap_depo = track.energy['depoTotal_no_ncap_depo']
        children = track.association['children']
        for childId in children:
            energy_no_ncap_depo += self.GetEnergyDepoWithDesendentsExcludeNCapEDepo(childId)
        return energy_no_ncap_depo

    #-------------------------
    def GetChargeDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        charge = track.energy['depoTotal_charge']
        charge_75keV = track.energy['depoTotal_charge_th_75keV']
        charge_500keV = track.energy['depoTotal_charge_th_500keV']
        children = track.association['children']
        for childId in children:
            charge += self.GetChargeDepoWithDesendents(childId)[0]
            charge_75keV += self.GetChargeDepoWithDesendents(childId)[1]
            charge_500keV += self.GetChargeDepoWithDesendents(childId)[2]
        return [charge, charge_75keV, charge_500keV]

    #-------------------------
    def GetChargeDepoWithDesendentsExcludeNCapEDepo(self, trkId):
        track = self.tracks[trkId]
        charge_no_ncap_depo = track.energy['depoTotal_charge_no_ncap_depo']
        charge_75keV_no_ncap_depo = track.energy['depoTotal_charge_th_75keV_no_ncap_depo']
        charge_500keV_no_ncap_depo = track.energy['depoTotal_charge_th_500keV_no_ncap_depo']
        children = track.association['children']
        for childId in children:
            charge_no_ncap_depo += self.GetChargeDepoWithDesendentsExcludeNCapEDepo(childId)[0]
            charge_75keV_no_ncap_depo += self.GetChargeDepoWithDesendentsExcludeNCapEDepo(childId)[1]
            charge_500keV_no_ncap_depo += self.GetChargeDepoWithDesendentsExcludeNCapEDepo(childId)[2]
        return [charge_no_ncap_depo, charge_75keV_no_ncap_depo, charge_500keV_no_ncap_depo]

    #-------------------------
    def GetChargeMBoxDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        charge = track.energy['depoTotal_charge_MBox']
        charge_75keV = track.energy['depoTotal_charge_MBox_th_75keV']
        charge_500keV = track.energy['depoTotal_charge_MBox_th_500keV']
        children = track.association['children']
        for childId in children:
            charge += self.GetChargeMBoxDepoWithDesendents(childId)[0]
            charge_75keV += self.GetChargeMBoxDepoWithDesendents(childId)[1]
            charge_500keV += self.GetChargeMBoxDepoWithDesendents(childId)[2]
        return [charge, charge_75keV, charge_500keV]

    #-------------------------
    def GetLightDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        light_avg_220PEpMeV = track.energy['depoTotal_light_avg_220PEpMeV']
        light_avg_180PEpMeV = track.energy['depoTotal_light_avg_180PEpMeV']
        light_avg_140PEpMeV = track.energy['depoTotal_light_avg_140PEpMeV']
        light_avg_100PEpMeV = track.energy['depoTotal_light_avg_100PEpMeV']
        light_avg_35PEpMeV  = track.energy['depoTotal_light_avg_35PEpMeV']
        children = track.association['children']
        for childId in children:
            light_avg_220PEpMeV += self.GetLightDepoWithDesendents(childId)[0]
            light_avg_180PEpMeV += self.GetLightDepoWithDesendents(childId)[1]
            light_avg_140PEpMeV += self.GetLightDepoWithDesendents(childId)[2]
            light_avg_100PEpMeV += self.GetLightDepoWithDesendents(childId)[3]
            light_avg_35PEpMeV  += self.GetLightDepoWithDesendents(childId)[4]
        return [light_avg_220PEpMeV, light_avg_180PEpMeV, light_avg_140PEpMeV, light_avg_100PEpMeV, light_avg_35PEpMeV]

    #-------------------------
    def GetLightDepoWithDesendentsExcludeNCapEDepo(self, trkId):
        track = self.tracks[trkId]
        light_avg_220PEpMeV_exclude_ncap_edep = track.energy['depoTotal_light_avg_220PEpMeV_exclude_ncap_edepo']
        light_avg_180PEpMeV_exclude_ncap_edep = track.energy['depoTotal_light_avg_180PEpMeV_exclude_ncap_edepo']
        light_avg_140PEpMeV_exclude_ncap_edep = track.energy['depoTotal_light_avg_140PEpMeV_exclude_ncap_edepo']
        light_avg_100PEpMeV_exclude_ncap_edep = track.energy['depoTotal_light_avg_100PEpMeV_exclude_ncap_edepo']
        light_avg_35PEpMeV_exclude_ncap_edep  = track.energy['depoTotal_light_avg_35PEpMeV_exclude_ncap_edepo']
        children = track.association['children']
        for childId in children:
            light_avg_220PEpMeV_exclude_ncap_edep += self.GetLightDepoWithDesendentsExcludeNCapEDepo(childId)[0]
            light_avg_180PEpMeV_exclude_ncap_edep += self.GetLightDepoWithDesendentsExcludeNCapEDepo(childId)[1]
            light_avg_140PEpMeV_exclude_ncap_edep += self.GetLightDepoWithDesendentsExcludeNCapEDepo(childId)[2]
            light_avg_100PEpMeV_exclude_ncap_edep += self.GetLightDepoWithDesendentsExcludeNCapEDepo(childId)[3]
            light_avg_35PEpMeV_exclude_ncap_edep  += self.GetLightDepoWithDesendentsExcludeNCapEDepo(childId)[4]
        return [light_avg_220PEpMeV_exclude_ncap_edep, light_avg_180PEpMeV_exclude_ncap_edep, light_avg_140PEpMeV_exclude_ncap_edep, light_avg_100PEpMeV_exclude_ncap_edep, light_avg_35PEpMeV_exclude_ncap_edep]

    #-------------------------
    def GetLightMBoxDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        light_avg_220PEpMeV = track.energy['depoTotal_light_avg_220PEpMeV_MBox']
        light_avg_180PEpMeV = track.energy['depoTotal_light_avg_180PEpMeV_MBox']
        light_avg_140PEpMeV = track.energy['depoTotal_light_avg_140PEpMeV_MBox']
        light_avg_100PEpMeV = track.energy['depoTotal_light_avg_100PEpMeV_MBox']
        light_avg_35PEpMeV  = track.energy['depoTotal_light_avg_35PEpMeV_MBox']
        children = track.association['children']
        for childId in children:
            light_avg_220PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[0]
            light_avg_180PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[1]
            light_avg_140PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[2]
            light_avg_100PEpMeV += self.GetLightMBoxDepoWithDesendents(childId)[3]
            light_avg_35PEpMeV  += self.GetLightMBoxDepoWithDesendents(childId)[4]
        return [light_avg_220PEpMeV, light_avg_180PEpMeV, light_avg_140PEpMeV, light_avg_100PEpMeV, light_avg_35PEpMeV]

    #-------------------------
    """def GetEnergyDepoWithAncestor(self, acId):
        energy = 0
        for track in self.tracks:
            ancestor = track.association['ancestor']
            if ancestor == acId:
                energy += track.energy['depoTotal']
        return energy"""

    #-----------------------
    def GetEnuFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-2]
        filename = filename.replace('MeV', '')
        return float(filename)  # Marley file names already in MeV

    #-----------------------
    def GetnuPDGFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-3]
        return filename

    #-----------------------
    def GetFileName(self):
        return self.simTree.GetFile().GetName()
# ------------------------
if __name__ == "__main__":
    event = Event(sys.argv[1])
    #event.Jump(0)
    #event.PrintVertex()
    #event.Next()
    #event.PrintVertex()
    #event.PrintTracks(0,8)
    # event.PrintTrack(1)
