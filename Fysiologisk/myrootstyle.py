from ROOT import gROOT,gPad,gStyle,TCanvas,TFile,TLine,TLatex,TAxis,TLegend,TPostScript
from ROOT import TH2D,TArrow,TCut,TPad,TPaveText,TGraph,TGraphAsymmErrors
from ROOT import TGraph2D,TStyle,TBranch,gSystem,gDirectory
from ROOT import TPave,TPaveStats


drStyle = TStyle("DR","DR style")

# use plain black on white colors
icol=0
drStyle.SetFrameBorderMode(icol)
drStyle.SetCanvasBorderMode(icol)
drStyle.SetPadBorderMode(icol)
drStyle.SetPadColor(icol)
drStyle.SetCanvasColor(icol)
drStyle.SetStatColor(icol)
#drStyle.SetFillColor(icol)

# set the paper & margin sizes
drStyle.SetPaperSize(20,26)
drStyle.SetPadTopMargin(0.05)
drStyle.SetPadRightMargin(0.05)
drStyle.SetPadBottomMargin(0.16)
drStyle.SetPadLeftMargin(0.12)

# use large fonts
font=42
tsize=0.05
drStyle.SetTextFont(font)


drStyle.SetTextSize(tsize)
drStyle.SetLabelFont(font,"x")
drStyle.SetTitleFont(font,"x")
drStyle.SetLabelFont(font,"y")
drStyle.SetTitleFont(font,"y")
drStyle.SetLabelFont(font,"z")
drStyle.SetTitleFont(font,"z")

drStyle.SetLabelSize(tsize,"x")
drStyle.SetTitleSize(tsize,"x")
drStyle.SetLabelSize(tsize,"y")
drStyle.SetTitleSize(tsize,"y")
drStyle.SetLabelSize(tsize,"z")
drStyle.SetTitleSize(tsize,"z")


#use bold lines and markers
drStyle.SetMarkerStyle(20)
drStyle.SetMarkerSize(1.2)
drStyle.SetHistLineWidth(2)
drStyle.SetLineStyleString(2,"[12 12]") # postscript dashes

#get rid of X error bars and y error bar caps
#drStyle.SetErrorX(0.001)

#do not display any of the standard histogram decorations
drStyle.SetOptTitle(0)
#drStyle.SetOptStat(1111)
drStyle.SetOptStat(0)
#drStyle.SetOptFit(1111)
drStyle.SetOptFit(0)

# put tick marks on top and RHS of plots
drStyle.SetPadTickX(1)
drStyle.SetPadTickY(1)

gROOT.SetStyle("Plain")

#gStyle.SetPadTickX(1)
#gStyle.SetPadTickY(1)
gROOT.SetStyle("DR")
gROOT.ForceStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
# overwrite dr styles

drStyle.SetMarkerSize(1.0)
drStyle.SetPadLeftMargin(0.14)
drStyle.SetPadRightMargin(0.03)
drStyle.SetPadBottomMargin(0.12)
drStyle.SetPadTopMargin(0.02)
drStyle.SetFrameFillColor(0)

def DRLabel(x,y,shift,Preliminary=False,color=1):
  l=TLatex()
  l.SetNDC()
  l.SetTextFont(72)
  l.SetTextColor(color)
  l.DrawLatex(x,y,"Danmarks Radio")
  if (Preliminary):
    p=TLatex()
    p.SetNDC()
    p.SetTextFont(42)
    p.SetTextColor(color)
    p.DrawLatex(x+shift,y,"Preliminary")



def DRVersion(version="1.0",x=0.88,y=0.975,color=1):
  if (version):
    l=TLatex()
    l.SetTextAlign(22)
    l.SetTextSize(0.04)
    l.SetNDC()
    l.SetTextFont(72)
    l.SetTextColor(color)
    l.DrawLatex(x,y,version)

def myText(x,y,color=1,size=0.08,text=""):
  l=TLatex()
  l.SetTextSize(size)
  l.SetNDC()
  l.SetTextColor(color)
  l.DrawLatex(x,y,text)


def DR_LABEL(x,y,color=1):
  l=TLatex()
  l.SetNDC()
  l.SetTextFont(72)
  l.SetTextColor(color)
  l.DrawLatex(x,y,"Danmarks Radio")
