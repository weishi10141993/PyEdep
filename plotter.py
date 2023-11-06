from event import Event
import sys
import matplotlib
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
        self.tt = np.array([])
        self.ee = np.array([])
        self.ll = np.array([])

        self.USER_COLORS = ['black', 'red', 'blue', 'magenta']

        nColor = len(self.USER_COLORS)        
        mm2m = 0.001
        mm2cm = 0.1
        for i, track in enumerate(self.event.tracks):   
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            for di in depoList:
                depo = self.event.depos[di]   
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                t = (depo.GetStart().T() + depo.GetStop().T()) / 2 # ns
                e = depo.GetEnergyDeposit() 
                l = depo.GetTrackLength() *mm2cm # most are 0.5 mm
                self.xx = np.append(self.xx, x)
                self.yy = np.append(self.yy, y)
                self.zz = np.append(self.zz, z)
                self.tt = np.append(self.tt, t)
                self.ee = np.append(self.ee, e)
                self.ll = np.append(self.ll, l)
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
    
    def Draw(self, axis='yz', value='time', markerSize=0.2, cmap='jet', vmax=2000):
        # particle, timing, dE/dx

        mapping = {'x': self.xx, 'y': self.yy, 'z': self.zz}
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5*2, 4), dpi=100)
        fig.suptitle(self.event.vertex.GetReaction())
        cb_ax = fig.add_axes([.94,.124,.03,.754])

        # particle plot
        ax1.scatter(mapping[axis[1]], mapping[axis[0]], c=self.cc, s=markerSize)
        if value == 'time':
            # timing plot
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.tt, cmap=cmap, vmin=1, vmax=vmax, norm=matplotlib.colors.LogNorm(), s=markerSize)
            cb_ax.set_xlabel('ns')

        elif value == 'charge':
            # charge plot
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ee*2, cmap=cmap, vmax=4, s=markerSize)
            cb_ax.set_xlabel('MeV/cm')

        
        fig.colorbar(plot_12, orientation='vertical', cax=cb_ax)
        # fig.tight_layout()

        for ax in (ax1, ax2):
            ax.set_ylim(-4, 4)
            ax.set_xlim(-2, 8)
            ax.tick_params(axis='y', direction='in', length=2)
            ax.tick_params(axis='x', direction='in', length=2)
            ax.set_xlabel(f'{axis[1]} [m]')
            if ax == ax1:
                ax.set_ylabel(f'{axis[0]} [m]')
        
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

            ax1.text(xpos, ypos, name, color=color)
            ypos -= 0.35

        plt.show()        

    def hist_dqdx(self):
        plt.hist(self.ee*2, range=(0,6), bins=100)
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

