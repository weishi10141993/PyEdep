from event import Event
import sys
import matplotlib.pyplot as plt
import numpy as np

class Plotter:

    def __init__(self, event):
        self.event = event
        self.Jump(0)
    
    def Collect(self):
        self.xx = np.array([])
        self.yy = np.array([])
        self.zz = np.array([])
        self.cc = np.array([])
        self.USER_COLORS = ['black', 'red', 'blue', 'magenta']

        nColor = len(self.USER_COLORS)        
        mm2m = 0.001
        for i, track in enumerate(self.event.tracks):   
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            for di in depoList:
                depo = self.event.depos[di]   
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                self.xx = np.append(self.xx, x)
                self.yy = np.append(self.yy, y)
                self.zz = np.append(self.zz, z)
                self.cc = np.append(self.cc, self.USER_COLORS[ancestor % nColor])

    def Jump(self, entryNo):
        self.event.Jump(entryNo)
        self.Collect()

    def Next(self):
        self.event.Next()
        self.Collect()
    
    def Prev(self):
        self.event.Prev()
        self.Collect()
    
    def Draw(self, axis='yz', markerSize=0.2):

        mapping = {'x': self.xx, 'y': self.yy, 'z': self.zz}
        
        fig, ax = plt.subplots()
        plt.scatter(mapping[axis[1]], mapping[axis[0]], c=self.cc, s=markerSize)
        fig.suptitle(self.event.vertex.GetReaction())
        plt.ylim(-4, 4)
        plt.xlim(-2, 8)
        plt.ylabel(f'{axis[0]} [m]')
        plt.xlabel(f'{axis[1]} [m]')

        ax.tick_params(axis='y', direction='in', length=2)
        ax.tick_params(axis='x', direction='in', length=2)

        xpos = -1.8
        ypos = 3.6
        nColor = len(self.USER_COLORS)        
        for i, particle in enumerate(self.event.vertex.Particles):   
            name = particle.GetName()
            color = self.USER_COLORS[i % nColor]
            # pdg = particle.GetPDGCode()
            # name = particle.GetName()
            # trkId = particle.GetTrackId()
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            name = '%s: %.1f MeV' % (name, KE)


            plt.text(xpos, ypos, name, color=color)
            ypos -= 0.35


        plt.show()        

    #---------------------------------------------
    def DrawROOT(self, dim2d='yz', markerSize=0.2):
        import ROOT
        from ROOT import TH2F, TMarker, TCanvas, TLatex
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetMarkerStyle(24)
        ROOT.gStyle.SetMarkerSize(markerSize)
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta]
        nColor = len(colors)

        c1 = TCanvas("c1", self.event.vertex.GetReaction(), 800, 800)

        # canvas of 5m x 5m
        dummy = TH2F("dummy", "", 100, -5, 5, 100, -5, 5)
        dummy.GetXaxis().SetTitle('[m]')
        dummy.GetYaxis().SetTitle('[m]')
        dummy.Draw()

        mm2m = 0.001
        txts = []
        markers = []
        txtX = 0.15
        txtY = 0.85
        txtSize = 0.03
        # nDepo = 0
        for i, track in enumerate(self.event.tracks):   
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            color = colors[ancestor % nColor]


            if ancestor == i:
                txt = TLatex()
                txt = txt.DrawLatexNDC(txtX, txtY, track.GetName())
                txt.SetTextColor(color)
                txt.SetTextSize(txtSize)
                txts.append(txt)
                txtY -= txtSize
                
            for j, di in enumerate(depoList):
                depo = self.event.depos[di]   
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                mapping = {'x': x, 'y': y, 'z': z}

                m = TMarker(mapping[dim2d[1]], mapping[dim2d[0]], 24)
                m.SetMarkerColor(color)

                m.Draw()
                markers.append(m)
                # nDepo += 1
        
        # print('depo points drawn: ', nDepo, '| total depo: ', self.event.depos.size)

        ROOT.gPad.Update()
        return c1

if __name__ == "__main__":
    event = Event(sys.argv[1])  
    p = Plotter(event)
    p.Next()
    p.Draw('yz')    
    # c1 = p.DrawROOT('xz', 0.2)
    # input('press a key to continue ...')        

