import ROOT
from ROOT import *
from array import array
import numpy
import math
from math import *
from collections import Counter


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inFile = TFile.Open("../out.root")
tree=inFile.Get("zlepTree")
nentries=tree.GetEntries()

stack=THStack("stack","stack")
stack.SetTitle("BR(Hinv)=10%, 900/fb, eLpR 80/30 polarization;recoil mass [GeV];Entries");
r0=TH1D("signal",";recoil mass [GeV];Entries",60,100.,160)
r1=TH1D("SMhiggs",";recoil mass [GeV];Entries",60,100.,160)
r2=TH1D("2fbkg",";recoil mass [GeV];Entries",60,100.,160)
r31=TH1D("31fbkg",";recoil mass [GeV];Entries",60,100.,160)
r32=TH1D("32fbkg",";recoil mass [GeV];Entries",60,100.,160)
r41=TH1D("41fbkg",";recoil mass [GeV];Entries",60,100.,160)
r42=TH1D("42fbkg",";recoil mass [GeV];Entries",60,100.,160)
r43=TH1D("43fbkg",";recoil mass [GeV];Entries",60,100.,160)
r44=TH1D("44fbkg",";recoil mass [GeV];Entries",60,100.,160)
r45=TH1D("45fbkg",";recoil mass [GeV];Entries",60,100.,160)

stackRL=THStack("stackRL","stackRL")
stackRL.SetTitle("BR(Hinv)=10%, 900/fb, eRpL 80/30 polarization;recoil mass [GeV];Entries");
r0RL=TH1D("signalRL",";recoil mass [GeV];Entries",60,100.,160)
r1RL=TH1D("SMhiggsRL",";recoil mass [GeV];Entries",60,100.,160)
r2RL=TH1D("2fbkgRL",";recoil mass [GeV];Entries",60,100.,160)
r31RL=TH1D("31fbkgRL",";recoil mass [GeV];Entries",60,100.,160)
r32RL=TH1D("32fbkgRL",";recoil mass [GeV];Entries",60,100.,160)
r41RL=TH1D("41fbkgRL",";recoil mass [GeV];Entries",60,100.,160)
r42RL=TH1D("42fbkgRL",";recoil mass [GeV];Entries",60,100.,160)
r43RL=TH1D("43fbkgRL",";recoil mass [GeV];Entries",60,100.,160)
r44RL=TH1D("44fbkgRL",";recoil mass [GeV];Entries",60,100.,160)
r45RL=TH1D("45fbkgRL",";recoil mass [GeV];Entries",60,100.,160)

for i in range(nentries) :
    tree.GetEntry(i)
    if i%10000==0:
	print "Event ",i
    if tree.nmbf ==0:
	if tree.iseLpR==1:
	    r0.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r0RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==1:
	if tree.iseLpR==1:
	    r1.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r1RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==2:
	if tree.iseLpR==1:
	    r2.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r2RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==31:
	if tree.iseLpR==1:
	    r31.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r31RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==32:
	if tree.iseLpR==1:
	    r32.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r32RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==41:
	if tree.iseLpR==1:
	    r41.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r41RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==42:
	if tree.iseLpR==1:
	    r42.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r42RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==43:
	if tree.iseLpR==1:
	    r43.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r43RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==44:
	if tree.iseLpR==1:
	    r44.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r44RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
    elif tree.nmbf==45:
	if tree.iseLpR==1:
	    r45.Fill(tree.m_rec,tree.eventWeight*tree.normalization)
	else:
	    r45RL.Fill(tree.m_rec,tree.eventWeight*tree.normalization)

stack.Add(r2)
stack.Add(r31)
stack.Add(r32)
stack.Add(r41)
stack.Add(r42)
stack.Add(r43)
stack.Add(r44)
stack.Add(r45)
stack.Add(r1)
stack.Add(r0)
stackRL.Add(r2RL)
stackRL.Add(r31RL)
stackRL.Add(r32RL)
stackRL.Add(r41RL)
stackRL.Add(r42RL)
stackRL.Add(r43RL)
stackRL.Add(r44RL)
stackRL.Add(r45RL)
stackRL.Add(r1RL)
stackRL.Add(r0RL)

r0.SetFillColor(1)
r1.SetFillColor(10)
r2.SetFillColor(2)
r31.SetFillColor(3)
r32.SetFillColor(4)
r41.SetFillColor(5)
r42.SetFillColor(6)
r43.SetFillColor(7)
r44.SetFillColor(8)
r45.SetFillColor(9)

r0RL.SetFillColor(1)
r1RL.SetFillColor(10)
r2RL.SetFillColor(2)
r31RL.SetFillColor(3)
r32RL.SetFillColor(4)
r41RL.SetFillColor(5)
r42RL.SetFillColor(6)
r43RL.SetFillColor(7)
r44RL.SetFillColor(8)
r45RL.SetFillColor(9)

print "eLpR signal integral: ",r0.Integral()
print "eLpR SM Higgs integral: ",r1.Integral()
print "eLpR 2f integral: ",r2.Integral()
print "eLpR 3f ap integral: ",r31.Integral()
print "eLpR 3f ea integral: ",r32.Integral()
print "eLpR 4f Wev integral: ",r41.Integral()
print "eLpR 4f WW integral: ",r42.Integral()
print "eLpR 4f Zee integral: ",r43.Integral()
print "eLpR 4f Zvv integral: ",r44.Integral()
print "eLpR 4f ZZ integral: ",r45.Integral()
print "eRpL signal integral: ",r0RL.Integral()
print "eRpL SM Higgs integral: ",r1RL.Integral()
print "eRpL 2f integral: ",r2RL.Integral()
print "eRpL 3f ap integral: ",r31RL.Integral()
print "eRpL 3f ea integral: ",r32RL.Integral()
print "eRpL 4f Wev integral: ",r41RL.Integral()
print "eRpL 4f WW integral: ",r42RL.Integral()
print "eRpL 4f Zee integral: ",r43RL.Integral()
print "eRpL 4f Zvv integral: ",r44RL.Integral()
print "eRpL 4f ZZ integral: ",r45RL.Integral()

leg = TLegend(0.7,0.7,0.9,0.9)
leg.SetHeader("eLpR","C")
leg.AddEntry(r0,"Hinv signal","f")
leg.AddEntry(r1,"SM Higgs","f")
leg.AddEntry(r2,"2f bkgs","f")
leg.AddEntry(r31,"3f ap","f")
leg.AddEntry(r32,"3f ea","f")
leg.AddEntry(r41,"4f Wev","f")
leg.AddEntry(r42,"4f WW","f")
leg.AddEntry(r43,"4f Zee","f")
leg.AddEntry(r44,"4f Zvv","f")
leg.AddEntry(r45,"4f ZZ","f")

canvas = TCanvas("can","can",600,600)
stack.Draw("hist")
leg.Draw()
canvas.SaveAs("m_rec_eLpR.png")

canvasRL = TCanvas("canRL","canRL",600,600)
stackRL.Draw("hist")
leg.SetHeader("eRpL","C")
leg.Draw()
canvas.SaveAs("m_rec_eRpL.png")
