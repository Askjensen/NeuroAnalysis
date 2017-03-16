# coding=utf-8
import os
import sys

import math
import warnings
from collections import namedtuple

import multiprocessing
import numpy as np
from numpy import (isscalar, r_, log, around, unique, asarray,
                   zeros, arange, sort, amin, amax, any, atleast_1d,
                   sqrt, ceil, floor, array, poly1d, compress,
                   pi, exp, ravel, angle, count_nonzero)
import pandas as pd
from array import array

import scipy
import scipy.stats as scistats
from ROOT import * #TH1F, TSpectrum, TFile, TCanvas, TLegend, gStyle, TAxis
#perhaps a transition to rootpy could be made entirely...:
#from rootpy.plotting import Hist

from bin.definitions import delimiter, folder_path, EventBinsPos, binscale, sigmapeaksinterval, peakamplitude, \
    EventBinsNames, f, dataset_event_names, filelist, Events_list, addconstant, test_length


def dataextract(files, label):
    fhandle = open(folder_path + files, 'r')

    temp_data = list()
    n = 0

    for p in fhandle:

        p1 = p.split(delimiter)

        if n == 0:
            label_lookup = p1.index(label)

        elif n != 0:

            if p1[label_lookup] != '':
                temp_data.append(float(p1[label_lookup].replace(',', '.')))

            elif p1[label_lookup] == '':
                if(len(temp_data)==0):
                    temp_data.append(0)
                else:
                    temp_data.append(temp_data[len(temp_data) - 1])

        n = n + 1

    fhandle.close()
    return (temp_data)


#create dir if not exist
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def createplotsFullRange(dataset_index_sub, dataarray,timewindow,name):
    timeinterval=12000 # total of 12000/60 = 200 minutes now - corrected later on
    #2D histograms:
    histphasicVSpeaksx = TH2F("histphasicVSpeaksx", "Phasic vs peak-position per respondent", (timeinterval + timewindow) / timewindow, 0,timeinterval,100,0,0.1)
    histtonicVSpeaksx = TH2F("histtonicVSpeaksx", "tonic vs peak-position per respondent", (timeinterval + timewindow) / timewindow, 0,timeinterval,100,0,0.1)
    histphasicVStonic = TH2F("histphasicVStonic", "Phasic vs tonic per respondent", 100,0,0.02,100,0,0.02)

    #stack
    hs = THStack("hs", "test stacked histograms")

    #1D hists
    histrawEDA = TH1F("histrawEDA", "raw EDA summed over nuuity normalised respondents",
                      (timeinterval + timewindow) / timewindow, 0, timeinterval)
    histtonicfull = TH1F("tonicsperminfull", "Tonic per respondent", (timeinterval + timewindow) / timewindow, 0,
                         timeinterval)
    histtonicfull.GetYaxis().SetTitle('Tonic component  / ' + str(timewindow) + ' sec.')
    histtonicfull.GetXaxis().SetTitle('Time [s]')
    histtonicfull.GetYaxis().SetTitleOffset(1.4)
    histtonicfull.SetStats(False)

    histphasicfull = TH1F("phasicsperminfull", "Phasic per respondent", (timeinterval + timewindow) / timewindow, 0,
                         timeinterval)
    histphasicfull.GetYaxis().SetTitle('Phasic component  / ' + str(timewindow) + ' sec.')
    histphasicfull.GetXaxis().SetTitle('Time [s]')
    histphasicfull.GetYaxis().SetTitleOffset(1.4)
    histphasicfull.SetStats(False)
    histpeaksfull = TH1F("histpeaksfull", "Peaks per respondent", (timeinterval + timewindow) / timewindow, 0,
                          timeinterval)
    histpeaksfull.GetYaxis().SetTitle('Peak component  / ' + str(timewindow) + ' sec.')
    histpeaksfull.GetXaxis().SetTitle('Time [s]')
    histpeaksfull.GetYaxis().SetTitleOffset(1.4)
    histpeaksfull.SetStats(False)
    histpeaksphasicamp = TH1F("peaksxampsperminfull", "Peaks*amplitude per respondent", (timeinterval + timewindow) / timewindow, 0,
                          timeinterval)
    histpeaksphasicamp.GetYaxis().SetTitle('Peak amplitude  / ' + str(timewindow) + ' sec.')
    histpeaksphasicamp.GetXaxis().SetTitle('Time [s]')
    histpeaksphasicamp.GetYaxis().SetTitleOffset(1.4)
    histpeaksphasicamp.SetStats(False)
    canvas = TCanvas("c", "Sum over respondents", 1200, 800)
    maxvalx = 0
    maxvaly = 0
    s = TSpectrum()
    #samples = zeros(len(dataset_index_sub))
    #iresp=0
    for respondent in dataset_index_sub:
        canvas1 = TCanvas("c1", "Count of peaks per respondent", 1200, 800)
        canvas1.cd()
        respdataarray = getrespondent(dataarray,respondent)

        hist = TH1F("hist" + respondent, "hist" + respondent,
                                    int(len(respdataarray.index)/20),
                                    respdataarray.index[0], respdataarray.index[-1])
        startpos = getfirstevent(Events_list, EventBinsNames, respondent)
        a = startpos + (1 / binscale)
        b = EventBinsPos[respondent][len(EventBinsPos[respondent]) - 1] - (1 / binscale)
        for i in range(0, len(respdataarray) - 1):
            hist.Fill(respdataarray.index[i], respdataarray.values[i])
            #if respdataarray.index[i] >= a and respdataarray.index[0] < b:
                #samples[iresp]+= respdataarray.values[i]
        hist.Sumw2()
        hist.GetXaxis().SetRangeUser(a, b)
        #hist.Scale(1./hist.Integral())
        hist.Draw()
        canvas1.Update()
        npeaks = s.Search(hist, sigmapeaksinterval, "noMarkov same nobackground", peakamplitude)
        bg = s.Background(hist, 20, "Compton same")
        canvas1.Update()
        for ibin in range(hist.GetNbinsX()):
            histrawEDA.SetBinContent(ibin,hist.GetBinContent(ibin)/hist.Integral()+histrawEDA.GetBinContent(ibin))

        #hist.Scale(1. / hist.Integral())
        # add constant to avoid negative values
        if addconstant:
            for i in range(0, hist.GetNbinsX()):
                hist.SetBinContent(i, hist.GetBinContent(i) + 100)
                bg.SetBinContent(i, bg.GetBinContent(i) + 100)
            orig_integral_of_full = hist.Integral()
            hist.Add(bg, -1)
            bg.Scale(1/orig_integral_of_full)
            hist.Scale(1. / (orig_integral_of_full-100))

        else:
            hist.Add(bg, -1)
            hist.Scale(1. / hist.Integral())
            bg.Scale(1. / bg.Integral())

        for x in range(0,npeaks):
            peakx = s.GetPositionX()[x]-a
            histpeaksfull.Fill(peakx * binscale)
            histpeaksphasicamp.Fill(peakx * binscale,hist.GetBinContent(hist.FindBin(peakx)))
            if((b-a)* binscale >maxvalx): maxvalx = (b - a) * binscale
            histphasicVSpeaksx.Fill(peakx * binscale, hist.GetBinContent(hist.FindBin(peakx)))
            histtonicVSpeaksx.Fill(peakx * binscale, bg.GetBinContent(bg.FindBin(peakx)))

        for ibin in range(bg.FindBin(a),bg.FindBin(b)):
            binx = bg.GetBinCenter(ibin)-a
            histtonicfull.Fill(binx * binscale, bg.GetBinContent(ibin))
            histphasicfull.Fill(binx * binscale, hist.GetBinContent(ibin))
            histphasicVStonic.Fill(bg.GetBinContent(ibin),hist.GetBinContent(ibin))
    if(histpeaksfull.GetMaximum()>maxvaly): maxvaly = histpeaksfull.GetMaximum()
    canvas.cd()

    print "correlation between phasic and number of peaks: " + str(histphasicVSpeaksx.GetCorrelationFactor())
    print "correlation between tonic and number of peaks: " + str(histtonicVSpeaksx.GetCorrelationFactor())
    print "correlation between phasic and tonic: " + str(histphasicVStonic.GetCorrelationFactor())
    print "covarience between phasic and number of peaks: " + str(histphasicVSpeaksx.GetCovariance())
    print "covarience between tonic and number of peaks: " + str(histtonicVSpeaksx.GetCovariance())
    print "covarience between phasic and tonic: " + str(histphasicVStonic.GetCovariance())
    #save2d
    histphasicVSpeaksx.GetYaxis().SetRangeUser(0, maxvaly*1.3)
    histphasicVSpeaksx.GetXaxis().SetRangeUser(0, maxvalx)
    histphasicVSpeaksx.Draw()
    canvas.Update()
    histphasicVSpeaksx.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/phasicVSpeaksx_' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/phasicVSpeaksx_' + name + 'overview' + '.png'))

    histtonicVSpeaksx.GetYaxis().SetRangeUser(0, maxvaly * 1.3)
    histtonicVSpeaksx.GetXaxis().SetRangeUser(0, maxvalx)
    histtonicVSpeaksx.Draw()
    canvas.Update()
    histtonicVSpeaksx.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/tonicVSpeaksx_' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/tonicVSpeaksx_' + name + 'overview' + '.png'))

    histphasicVStonic.GetYaxis().SetRangeUser(0, maxvaly * 1.3)
    histphasicVStonic.GetXaxis().SetRangeUser(0, maxvalx)
    histphasicVStonic.Draw()
    canvas.Update()
    histphasicVStonic.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/phasicVStonic_' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/phasicVStonic_' + name + 'overview' + '.png'))

    #save 1D
    histpeaksfull.GetYaxis().SetRangeUser(0, maxvaly*1.3)
    histpeaksfull.GetXaxis().SetRangeUser(0, maxvalx)
    histpeaksfull.Draw()
    #leg.Draw()
    canvas.Update()
    histpeaksfull.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/timedistributionof' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/timedistributionof_' + name + 'overview' + '.png'))

    histrawEDA.GetYaxis().SetRangeUser(0, histpeaksphasicamp.GetMaximum()*1.3)
    histrawEDA.GetXaxis().SetRangeUser(0, maxvalx)
    histrawEDA.Draw()
    canvas.Update()
    histrawEDA.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/timedistributionofrawEDAsignal' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/timedistributionofrawEDAsignal' + name + 'overview' + '.png'))

    histpeaksphasicamp.GetYaxis().SetRangeUser(0, histpeaksphasicamp.GetMaximum()*1.3)
    histpeaksphasicamp.GetXaxis().SetRangeUser(0, maxvalx)
    histpeaksphasicamp.Draw()
    #leg.Draw()
    canvas.Update()
    histpeaksphasicamp.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/timedistributionofpeaksxphasic' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/ctimedistributionofpeaksxphasic' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/timedistributionofpeaksxphasic' + name + 'overview' + '.png'))

    histphasicfull.GetXaxis().SetRangeUser(0, maxvalx)
    maxvaly=histphasicfull.GetMaximum()
    histphasicfull.GetYaxis().SetRangeUser(0, maxvaly*1.3)
    histphasicfull.Draw()
    #leg.Draw()
    canvas.Update()
    histphasicfull.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/timedistributionofphasic' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/ctimedistributionof_phasic' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/timedistributionof_phasic' + name + 'overview' + '.png'))

    maxvaly=histtonicfull.GetMaximum()
    histtonicfull.GetXaxis().SetRangeUser(0, maxvalx)
    histtonicfull.GetYaxis().SetRangeUser(0, maxvaly*1.3)
    histtonicfull.Draw()
    #leg.Draw()
    canvas.Update()
    histtonicfull.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/timedistributionoftonic' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/ctimedistributionof_tonic' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/timedistributionof_tonic' + name + 'overview' + '.png'))

    histtonicfull.SetFillColor(1)
    histphasicfull.SetFillColor(2)
    histtonicfull.SetFillStyle(3001)
    histphasicfull.SetFillStyle(3001)
    hs.Add(histtonicfull)
    hs.Add(histphasicfull)
    hs.Draw()
    canvas.Update()
    canvas.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/stacked_' + name + 'overview' + '.root'))
    canvas.SaveAs(os.path.join(sys.path[0],'../out/results/stacked_' + name + 'overview' + '.png'))
    canvas.Close()



def phasic_component(dataset_index_sub, dataarray, window):
    values = list()

    for p in dataset_index_sub:
        tonic_temp = (pd.rolling_mean(dataarray.loc[p], window))

        returndata1 = abs(((pd.Series.subtract(tonic_temp, dataarray.loc[p])) * (
            -1.0)))  # ABS: Negative peaks transformed into positive peaks

        values.extend(returndata1.tolist())

    returndata = pd.Series(values, index=dataarray.index)

    return (returndata)


def tonic_component(dataset_index_sub, dataarray, window):
    values = list()

    for p in dataset_index_sub:
        tonic_temp = (pd.rolling_mean(dataarray.loc[p], window))
        values.extend(tonic_temp)

    return pd.Series(values, index=dataarray.index)


def getrespondent(dataarray,respondent):
    try:
        respdataarray = dataarray.loc[respondent]
    except Exception, e:
        print respondent + "was not in dataarray index - skipping respondent for plots"
        print >> f, respondent + "was not in dataarray index - skipping respondent for plots - error: " + e
        return
    return respdataarray


def normalize_series(dataset_index_sub, dataseries):
    list_of_dataframes = []

    for p in dataset_index_sub:
        dataseries1 = (dataseries.loc[p] - dataseries.loc[p].min()) / (dataseries.loc[p].max() - dataseries.loc[p].min())
        dataseries1['Names'] = pd.Series(p, index=dataseries1.index)
        dataseries1['Syncpos']=dataseries1.index
        dataseries1.set_index(['Names', 'Syncpos'], inplace=True)

        list_of_dataframes.append(dataseries1)
    out = pd.concat(list_of_dataframes)
    return out


def findpeaksandstore(canvas, comparisonlist, dataarray, dataset_index_sub, ievent, ikey, leg, maxvalx, maxvaly, peakarray,
                      phasicarray, tonicarray):
    for event in comparisonlist[ikey]:
        for respondent in dataset_index_sub:
            respdataarray = getrespondent(dataarray, respondent)
            try:
                respdataarray
            except NameError:
                continue

            for clipindex in range(0, len(EventBinsNames[respondent]) - 1):
                if EventBinsNames[respondent][clipindex] != event: continue
                if (clipindex % 2 == 0):
                    canvas1 = TCanvas("c1", "Count of peaks per respondent", 1200, 800)
                    canvas1.cd()
                    s = TSpectrum()
                    hist = TH1F("hist" + respondent + event, "hist" + respondent + event,
                                int(len(respdataarray.index) / 20),
                                respdataarray.index[0], respdataarray.index[-1])
                    a = EventBinsPos[respondent][clipindex] + (1 / binscale)
                    b = EventBinsPos[respondent][clipindex + 1] - (1 / binscale)
                    nbins = (hist.FindBin(b) - hist.FindBin(a))
                    # phasic = TH1F("histphasic" + respondent + event, "hist" + respondent + event,0,b-a, nbins)
                    for i in range(0, len(respdataarray) - 1):
                        hist.Fill(respdataarray.index[i], respdataarray.values[i])
                        # if respdataarray.index[i]>= a:
                        # if respdataarray.index[i]<=b:
                        # phasic.Fill(respdataarray.index[i], respdataarray.values[i])
                    hist.GetXaxis().SetRangeUser(a, b)
                    if (hist.Integral(hist.FindBin(a), hist.FindBin(b) != 0)):
                        hist.Scale(1. / hist.Integral(hist.FindBin(a), hist.FindBin(b)))
                    # phasic.GetXaxis().SetRangeUser(a, b)
                    hist.Draw()
                    canvas1.Update()
                    npeaks = s.Search(hist, sigmapeaksinterval, "noMarkov same nobackground", peakamplitude)
                    bg = s.Background(hist, 20, "Compton same")
                    canvas1.Update()
                    hist.Add(bg, -1)
                    for x in range(0, npeaks):
                        peakx = s.GetPositionX()[x] - a
                        peakarray[ievent].Fill(peakx * binscale)
                        if ((b - a) * binscale > maxvalx): maxvalx = (b - a) * binscale
                    for ibin in range(bg.FindBin(a), bg.FindBin(b)):
                        binx = bg.GetBinCenter(ibin) - a
                        tonicarray[ievent].Fill(binx * binscale, bg.GetBinContent(ibin))
                        phasicarray[ievent].Fill(binx * binscale, hist.GetBinContent(ibin))
        canvas.cd()
        peakarray[ievent].SetLineColor(ievent + 1)
        peakarray[ievent].SetFillColor(ievent + 1)
        peakarray[ievent].SetFillStyle(3003 + ievent)
        tonicarray[ievent].SetLineColor(ievent + 1)
        tonicarray[ievent].SetFillColor(ievent + 1)
        tonicarray[ievent].SetFillStyle(3003 + ievent)
        phasicarray[ievent].SetLineColor(ievent + 1)
        phasicarray[ievent].SetFillColor(ievent + 1)
        phasicarray[ievent].SetFillStyle(3003 + ievent)

        leg.AddEntry(peakarray[ievent], event, "l")
        canvas.Update()
        if (peakarray[ievent].GetMaximum() > maxvaly): maxvaly = peakarray[ievent].GetMaximum()
        ievent += 1
    return maxvalx, maxvaly


def npeaksspecificminutes(comparisonlist, dataarray, dataset_index_sub, name, timewindow):
    timeinterval = 12000  # total of 12000/60 = 200 minutes now - corrected later on
    histtonicfull = TH1F("tonicsperminfull", "Tonic Component", (timeinterval + timewindow) / timewindow, 0,
                         timeinterval)
    histtonicfull.GetYaxis().SetTitle('Arousal  / ' + str(timewindow) + ' sec.')
    histtonicfull.GetXaxis().SetTitle('Time [s]')
    histtonicfull.GetYaxis().SetTitleOffset(1.4)
    histtonicfull.SetStats(False)
    histphasicfull = TH1F("phasicperminfull", "Phasic component", (timeinterval + timewindow) / timewindow, 0,
                          timeinterval)
    histphasicfull.GetYaxis().SetTitle('Arousal  / ' + str(timewindow) + ' sec.')
    histphasicfull.GetXaxis().SetTitle('Time [s]')
    histphasicfull.GetYaxis().SetTitleOffset(1.4)
    histphasicfull.SetStats(False)
    histpeaksfull = TH1F("peaksperminfull", "Peaks in arousel", (timeinterval + timewindow) / timewindow, 0,
                         timeinterval)
    histpeaksfull.GetYaxis().SetTitle('Number of peaks  / ' + str(timewindow) + ' sec.')
    histpeaksfull.GetXaxis().SetTitle('Time [s]')
    histpeaksfull.GetYaxis().SetTitleOffset(1.4)
    histpeaksfull.SetStats(False)
    for ikey in comparisonlist.keys():
        canvas = TCanvas("c", "Count of peaks per respondent", 1200, 800)
        canvas.cd()
        leg = TLegend(0.5, 0.67, 0.88, 0.88)
        ievent = 0
        histpeaks = TH1F("peakspermin", "peaks per respondent", (timeinterval + timewindow) / timewindow, 0,
                         timeinterval)
        histpeaks.GetYaxis().SetTitle('Number of peaks  / ' + str(timewindow) + ' sec.')
        histpeaks.GetXaxis().SetTitle('Time [s]')
        histpeaks.GetYaxis().SetTitleOffset(1.4)
        histpeaks.SetStats(False)
        histphasic = TH1F("phasicspermin", "Phasic per respondent", (timeinterval + timewindow) / timewindow, 0,
                          timeinterval)
        histphasic.GetYaxis().SetTitle('Phasic component  / ' + str(timewindow) + ' sec.')
        histphasic.GetXaxis().SetTitle('Time [s]')
        histphasic.GetYaxis().SetTitleOffset(1.4)
        histphasic.SetStats(False)
        histtonic = TH1F("tonicspermin", "Tonic per respondent", (timeinterval + timewindow) / timewindow, 0,
                         timeinterval)
        histtonic.GetYaxis().SetTitle('Tonic component  / ' + str(timewindow) + ' sec.')
        histtonic.GetXaxis().SetTitle('Time [s]')
        histtonic.GetYaxis().SetTitleOffset(1.4)
        histtonic.SetStats(False)
        peakarray = []
        phasicarray = []
        tonicarray = []
        maxvalx = 0
        maxvaly = 0
        peakappend(comparisonlist, ikey, peakarray, phasicarray, timeinterval, timewindow, tonicarray)

        maxvalx, maxvaly = findpeaksandstore(canvas, comparisonlist, dataarray, dataset_index_sub, ievent, ikey, leg, maxvalx,
                                             maxvaly, peakarray, phasicarray, tonicarray)
        canvas.cd()
        histpeaks.GetYaxis().SetRangeUser(0, maxvaly * 1.3)
        histpeaks.GetXaxis().SetRangeUser(0, maxvalx)
        histpeaks.Draw()
        for ievent in range(len(comparisonlist[ikey]) - 1, -1, -1):
            peakarray[ievent].Draw("same")
        leg.Draw()
        canvas.Update()
        canvas.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/timedistributionof' + name + 'inseq' + ikey + '.root'))
        canvas.SaveAs(os.path.join(sys.path[0],'../out/results/timedistributionof' + name + 'inseq' + ikey + '.png'))

        histphasic.GetXaxis().SetRangeUser(0, maxvalx)
        histphasic.Draw()
        maxvaly = 0
        for ievent in range(len(comparisonlist[ikey]) - 1, -1, -1):
            if maxvaly <= phasicarray[ievent].GetMaximum():
                maxvaly = phasicarray[ievent].GetMaximum()
            phasicarray[ievent].Draw("same")
        histphasic.GetYaxis().SetRangeUser(0, maxvaly * 1.3)
        leg.Draw()
        canvas.Update()
        canvas.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/timedistributionofphasic' + name + 'inseq' + ikey + '.root'))
        canvas.SaveAs(os.path.join(sys.path[0],'../out/results/timedistributionofphasic' + name + 'inseq' + ikey + '.png'))

        maxvaly = 0
        histtonic.GetXaxis().SetRangeUser(0, maxvalx)
        histtonic.Draw()
        for ievent in range(len(comparisonlist[ikey]) - 1, -1, -1):
            if maxvaly <= tonicarray[ievent].GetMaximum():
                maxvaly = tonicarray[ievent].GetMaximum()
            tonicarray[ievent].Draw("same")
        histtonic.GetYaxis().SetRangeUser(0, maxvaly * 1.3)
        leg.Draw()
        canvas.Update()
        canvas.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/timedistributionoftonic' + name + 'inseq' + ikey + '.root'))
        canvas.SaveAs(os.path.join(sys.path[0],'../out/results/timedistributionoftonic' + name + 'inseq' + ikey + '.png'))
        # print >> f, "\n KS probability including normalisation of " + str(comparisonlist[ikey][0]) + " and " + str(comparisonlist[ikey][1])
        # print >> f, histarray[0].KolmogorovTest(histarray[1], "OUN")
        # print >> f, "\n KS Comparison of " + str(comparisonlist[ikey][0]) + " and " + str(comparisonlist[ikey][1])
        # print >> f, histarray[0].KolmogorovTest(histarray[1], "OU")
        # print >> f,"\n Chi2 probability Comparison of " + str(comparisonlist[ikey][0]) + " and " + str(comparisonlist[ikey][1])
        # print >> f, histarray[0].Chi2Test(histarray[1],"UU")
        # print >> f, "\n mean " + str(comparisonlist[ikey][0])+": " + str(histarray[0].GetMean())
        # print >> f, "mean " + str(comparisonlist[ikey][1])+": " + str(histarray[0].GetMean())
        canvas.Close()

def peakappend(comparisonlist, ikey, peakarray, phasicarray, timeinterval, timewindow, tonicarray):
    for event in comparisonlist[ikey]:
        tmphist = TH1F("peakspermin" + event, event, (timeinterval + timewindow) / timewindow, 0, timeinterval)
        tmphist.GetXaxis().SetTitle('Time [s]')
        tmphist.GetYaxis().SetTitle('Number of peaks / ' + str(timewindow) + ' sec.')
        tmphist.SetStats(False)
        peakarray.append(tmphist)
        ph = TH1F("ph" + event, event, (timeinterval + timewindow) / timewindow, 0, timeinterval)
        phasicarray.append(ph)
        th = TH1F("th" + event, event, (timeinterval + timewindow) / timewindow, 0, timeinterval)
        tonicarray.append(th)
        return

def loop(canvas, comparisonlist, dataarray, dataset_index_sub, histarray, ievent, ikey, leg, maxval, name):
    for event in comparisonlist[ikey]:
        tmphist = TH1F("peaksperres" + event, event, 50, 0, 50)
        tmphist.GetXaxis().SetTitle('Number of peaks in sequence')
        tmphist.GetYaxis().SetTitle('Count')
        tmphist.SetStats(False)
        tmphist.SetFillStyle(3002)
        histarray.append(tmphist)
    for event in comparisonlist[ikey]:
        for respondent in dataset_index_sub:
            respdataarray = getrespondent(dataarray, respondent)

            hist = TH1F(str(respondent) + "hist", str(respondent) + "hist",
                        int((respdataarray.index[-1] - respdataarray.index[0]) * binscale),
                        respdataarray.index[0], respdataarray.index[-1])
            for i in range(0, len(respdataarray) - 1):
                hist.Fill(respdataarray.index[i], respdataarray.values[i])
            for clipindex in range(0, len(EventBinsPos[respondent]) - 1):
                if EventBinsNames[respondent][clipindex] != event: continue
                if (clipindex % 2 == 0):
                    s = TSpectrum()
                    # canvas1 = TCanvas("c1", "Count of peaks per respondent", 1200, 800)
                    # canvas1.cd()
                    hist.GetXaxis().SetRangeUser(EventBinsPos[respondent][clipindex],
                                                 EventBinsPos[respondent][clipindex + 1])
                    # hist.Draw()
                    # canvas1.Update()
                    np = s.Search(hist, sigmapeaksinterval, "noMarkov same nobackground", peakamplitude)
                    histarray[ievent].Fill(np)
                    # canvas1.Delete()
        canvas.cd()
        histarray[ievent].Draw("same")
        histarray[ievent].SetLineColor(ievent + 1)
        histarray[ievent].SetFillColor(ievent + 1)
        histarray[ievent].SetFillStyle(3003 + ievent)
        leg.AddEntry(histarray[ievent], event, "l")
        canvas.Update()
        if (histarray[ievent].GetMaximum() > maxval): maxval = histarray[ievent].GetMaximum()
        ievent += 1
    return maxval


def npeaksspecific(comparisonlist, dataarray, dataset_index_sub, name):
    for ikey in comparisonlist.keys():
        canvas = TCanvas("c", "Count of peaks per respondent", 1200, 800)
        canvas.cd()
        leg = TLegend(0.5, 0.67, 0.88, 0.88)
        ievent = 0
        histogram = TH1F("peaksperres", "peaks per respondent", 50, 0, 50)
        histogram.GetXaxis().SetTitle('Number of peaks in sequence')
        histogram.GetYaxis().SetTitle('Count')
        histogram.GetYaxis().SetTitleOffset(1.4)
        histogram.SetStats(False)
        histogram.Draw()
        canvas.Update()
        histarray = []
        maxval = 0
        maxval = loop(canvas, comparisonlist, dataarray, dataset_index_sub, histarray, ievent, ikey, leg, maxval, name)
        canvas.cd()

        histogram.GetYaxis().SetRangeUser(0, maxval)
        canvas.Update()
        # histogram.Write()
        leg.Draw()
        # canvas.Draw()
        # print >> f, "\n KS probability including normalisation of " + str(comparisonlist[ikey][0]) + " and " + str(
        #     comparisonlist[ikey][1])
        # print >> f, histarray[0].KolmogorovTest(histarray[1], "OUN")
        # print >> f, "\n KS Comparison of " + str(comparisonlist[ikey][0]) + " and " + str(comparisonlist[ikey][1])
        # print >> f, histarray[0].KolmogorovTest(histarray[1], "OU")
        # print >> f, "\n Chi2 probability Comparison of " + str(comparisonlist[ikey][0]) + " and " + str(
        #     comparisonlist[ikey][1])
        # print >> f, histarray[0].Chi2Test(histarray[1], "UU")
        # print >> f, "\n mean " + str(comparisonlist[ikey][0]) + ": " + str(histarray[0].GetMean()) + "standard deviation: " + str(histarray[0].GetRMS())
        # print >> f, "mean " + str(comparisonlist[ikey][1]) + ": " + str(histarray[1].GetMean()) + "standard deviation: " + str(histarray[1].GetRMS())
        #canvas.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/numberof' + name + 'inseq' + ikey + '.root'))
        histogram.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/numberof' + name + 'inseq' + ikey + '.root'))
        canvas.SaveAs(os.path.join(sys.path[0],'../out/results/numberof' + name + 'inseq' + ikey + '.png'))
        histogram.Delete()
        canvas.Close()


def meaninterval(dataset_index_sub, dataarray, comparisonlist,timewindow,name):
    timeinterval=12000 # total of 12000/60 = 200 minutes now - corrected later on
    for ikey in comparisonlist.keys():
        canvas = TCanvas("c", "mean EDA per respondent", 1200, 800)
        canvas.cd()
        leg = TLegend(0.5, 0.67, 0.88, 0.88)
        ievent = 0
        histogram = TH1F("meanedapertime", "mean EDA per respondent", (timeinterval+timewindow)/timewindow, 0, timeinterval)
        histogram.GetYaxis().SetTitle('Mean EDA of ' + name + 'data  / '+ str(timewindow) + ' sec.')
        histogram.GetXaxis().SetTitle('Time [s]')
        histogram.GetYaxis().SetTitleOffset(1.4)
        histogram.SetStats(False)
        histogram.Draw()
        canvas.Update()
        histarray = []
        maxvalx = 0
        maxvaly = 0
        for event in comparisonlist[ikey]:
            tmphist = TH1F("peakspermin" + event, event, (timeinterval+timewindow)/timewindow, 0, timeinterval)
            tmphist.SetStats(False)
            tmphist.SetFillStyle(3002)
            histarray.append(tmphist)

        for event in comparisonlist[ikey]:
            for respondent in dataset_index_sub:
                respdataarray = getrespondent(dataarray, respondent)

                for clipindex in range(0, len(EventBinsNames[respondent]) - 1):
                    if EventBinsNames[respondent][clipindex] != event: continue
                    if (clipindex % 2 == 0):
                        a = EventBinsPos[respondent][clipindex] * binscale
                        b = EventBinsPos[respondent][clipindex + 1] * binscale
                        for interval in range(0, int(round( (b-a)/ timewindow))):
                            hist = TH1F(str(respondent) + str(clipindex) + "hist", str(respondent) + "hist",
                                        int(len(respdataarray.values)), respdataarray.values.min(),
                                        respdataarray.values.max())
                            for i in range(0, len(respdataarray.values)):
                                if respdataarray.index[i] * binscale > a + interval * timewindow:
                                    if respdataarray.index[i] * binscale <= a + interval * timewindow + timewindow:
                                        hist.Fill(respdataarray.values[i])
                            meaninterval = hist.GetMean()
                            histarray[ievent].Fill(interval * timewindow, meaninterval)
                            if (interval * timewindow) > maxvalx: maxvalx = interval * timewindow
            histarray[ievent].SetLineColor(ievent+1)
            histarray[ievent].SetFillColor(ievent+1)
            histarray[ievent].SetFillStyle(3003+ievent)
            leg.AddEntry(histarray[ievent], event, "l")
            canvas.Update()
            if(histarray[ievent].GetMaximum()>maxvaly): maxvaly = histarray[ievent].GetMaximum()
            ievent += 1
        canvas.cd()
        histogram.GetYaxis().SetRangeUser(0, maxvaly*1.3)
        histogram.GetXaxis().SetRangeUser(0, maxvalx)
        histogram.Draw()
        for ievent in range(0,len(comparisonlist[ikey])):
            histarray[ievent].Draw("same")
        leg.Draw()
        canvas.Update()
        #histogram.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/compinterval'+name+ikey+'.root'))
        canvas.SaveAs(os.path.join(sys.path[0],'../out/results/compinterval'+name+ikey+'.root'))
        canvas.Close()
        return


def testsettings(dataset_index_sub, dataarray):
    threspeaks = TH2F("threspeaks", "peaks in baseline as a function og #sigma_{thres} and thres_{amplitude}", 5, .5, 5.5, 5, 0.05,0.3)
    threspeaksnomarkov = TH2F("threspeaksnomarkov", "peaks in baseline as a function og #sigma_{thres} and thres_{amplitude} with noMarkov", 5, .5, 5.5, 5, 0.05,0.3)
    peakposition = TH1F("peakpos", "peak position",36,0,180)
    peakstimulidiff = TH2F("peakdiff", "x_{peak} - x_{stimuli} (per respondent) as a function og #sigma_{thres} and thres_{amplitude}",5, .5, 5.5, 5, 0.05,0.3)
    npeakstimulidiff = TH2F("npeakstimulidiff", "#sum(x_{peak} - x_{stimuli})*(npeaks_{tot}-npeaks_{expected}) (per respondent) as a function og #sigma_{thres} and thres_{amplitude}",5, .5, 5.5, 5, 0.05,0.3)
    npeaks2d = TH2F("npeaks2d", "number of peaks (per respondent) as a function og #sigma_{thres} and thres_{amplitude}",5, .5, 5.5, 5, 0.05,0.3)
    npeaks1d = TH1F("npeaks1d", "number of peaks per respondent",21,0,20)
    threspeaks.GetXaxis().SetTitle("#sigma_{thres}")
    threspeaksnomarkov.GetXaxis().SetTitle("#sigma_{thres}")
    threspeaks.GetYaxis().SetTitle("thres_{amplitude}")
    threspeaksnomarkov.GetYaxis().SetTitle("thres_{amplitude}")
    threspeaks.GetZaxis().SetTitle("Peaks per perticipant")
    threspeaksnomarkov.GetZaxis().SetTitle("Peaks per perticipant")

    nres = len(dataset_index_sub)
    for respondent in dataset_index_sub:
        respdataarray = getrespondent(dataarray, respondent)
        stimulipositions = [0,14,25,35,48,60,70,80,95,106]

        hist = TH1F(str(respondent) + "hist", str(respondent) + "hist", int((respdataarray.index[-1]*binscale)),
                    respdataarray.index[0], respdataarray.index[-1])

        for i in range(0, len(respdataarray) - 1):
            hist.Fill(respdataarray.index[i], respdataarray.values[i])
        #orig_integral = hist.Integral()
        # loop over series in timeunit of size: peaks_window
        #hist.GetXaxis().SetTitle('#sigma_{thres}')
        #hist.GetYaxis().SetTitle('number of peaks')
        first = True
        for sigma in range(1,6):
            for thres in [0.051,0.101,0.151,0.201,0.251]:
                s = TSpectrum()
                np = s.Search(hist, sigma, "nodraw nobackground", thres)
                threspeaks.Fill(sigma,thres,np/float(nres))
                np = s.Search(hist, sigma, "noMarkov nodraw nobackground", thres)
                npdev = math.fabs(9 - np) + 1
                threspeaksnomarkov.Fill(sigma, thres, np / float(nres))
                npeaks2d.Fill(sigma, thres, float(np)/float(nres))
                npeaks1d.Fill(float(np),1/float(nres))
                for x in range(0, np):
                    peakx = s.GetPositionX()[x]
                    peakposition.Fill(peakx * binscale,1/float(nres))
                    peakstimulidiff.Fill(sigma, thres, getmindiff(peakx * binscale,stimulipositions)/float(nres))
                    npeakstimulidiff.Fill(sigma, thres, getmindiff(peakx * binscale,stimulipositions)*float(npdev)/float(nres)) #((1+float(np)-10)*float(nres)))

    threspeaks.Write()
    threspeaksnomarkov.Write()
    peakposition.Write()
    peakstimulidiff.Write()
#    peakstimulidiff.SaveAs('out/results/peakstimulidiff.png')
    npeaks2d.Write()
    npeaks1d.Write()
    npeakstimulidiff.Write()

def getmindiff(x,xstimuli):
    minimum = 100000
    for i in xstimuli:
        if math.fabs(x-i) < minimum:
            minimum = x-i
    return minimum

def npeaks(dataset_index_sub, dataarray,name,type):
    nbins = (len(dataset_event_names[0]) / 2)
    peakshist = TH1F("peaksperseq", "peaks per sequence", nbins, 0, nbins)
    peakshist.Sumw2()
    ibin = 0
    labels = sorted(dataset_event_names[0])
    for bin in range(0, nbins * 2):
        if bin % 2 == 0:
            ibin += 1
            peakshist.GetXaxis().SetBinLabel(ibin, labels[bin])

    for respondent in dataset_index_sub:
        respdataarray = getrespondent(dataarray, respondent)

        hist = TH1F(str(respondent) + "hist", str(respondent) + "hist", int((respdataarray.index[-1]*binscale)),
                    respdataarray.index[0], respdataarray.index[-1])

        for i in range(0, len(respdataarray) - 1):
            hist.Fill(respdataarray.index[i], respdataarray.values[i])
        orig_integral = hist.Integral()
        # loop over series in timeunit of size: peaks_window
        #hist.GetXaxis().SetTitle('Time-position [ms]')
        #hist.GetYaxis().SetTitle('EDA')
        # for sigma in range(1,5):
        #     for thres in [0.05]:  # [0.01,0.05,0.10]:
        #         for res in [1]:  # range(1,5)
        #             s = TSpectrum()
        #             s.SetResolution(res)
        #             cvs = TCanvas("c", "c", 1200, 800)
        #             cvs.cd()
        #             hist.Draw()
        #             cvs.Update()
        #             hist.Write()
        #             s.Search(hist, sigma, "noMarkov nodraw", thres)
        #             cvs.Update()
        #             bg = s.Background(hist, 20, "Compton nodraw")
        #             currentsettings = respondent + 'sigma' + str(sigma) + 'thres' + str(thres) + 'res' + str(res)
        #             cvs.Update()
        #             cvs.SaveAs(os.path.join(os.getcwd(), './../out/peaks/noMarkov' + currentsettings + '.png'))
        #             cvs.Close()
        #
        #             cmark = TCanvas("c", "c", 1200, 800)
        #             cmark.cd()
        #             hist.Draw()
        #             s.Search(hist, sigma, "nobackground nodraw", thres)
        #             cmark.Update()
        #             bg = s.Background(hist, 20, "nodraw nosmoothing")
        #             cmark.Update()
        #             cmark.SaveAs(os.path.join(os.getcwd(), './../out/peaks/nobackground_' + currentsettings + '.png'))
        #             cmark.Close()

        ## First define peaks for full sequence from baseline untill last sequence
        startpos = getfirstevent(Events_list,EventBinsNames,respondent)
        hist.GetXaxis().SetRangeUser(startpos + (1 / binscale),
                                     EventBinsPos[respondent][len(EventBinsPos[respondent])-1] - (1 / binscale))
        # hist.Scale(1./hist.Integral(hist.GetXaxis().FindBin(EventBinsPos[respondent][clipindex]),
        # hist.GetXaxis().FindBin(EventBinsPos[respondent][clipindex + 1])))
        cmark = TCanvas("c", "c", 1200, 800)
        cmark.cd()
        hist.Draw()
        s = TSpectrum()
        np = s.Search(hist, sigmapeaksinterval, "noMarkov same nobackground", peakamplitude)
        bg = s.Background(hist, 20, "Compton same")
        cmark.Update()
        dir = os.path.join(os.getcwd(), './../out/peaks' + name + '/overview/')
        ensure_dir(dir)
        cmark.SaveAs(dir + 'OverviewOf'+name+'Peaks_' + respondent + '.png')
        cmark.Close()

        ## then for individual "clips"
        for clipindex in range(0, len(EventBinsPos[respondent]) - 1):
            if (clipindex % 2 == 0):

                hist.GetXaxis().SetRangeUser(EventBinsPos[respondent][clipindex]+(1/binscale),
                                         EventBinsPos[respondent][clipindex + 1]-(1/binscale))
                #hist.Scale(1./hist.Integral(hist.GetXaxis().FindBin(EventBinsPos[respondent][clipindex]),
                                            #hist.GetXaxis().FindBin(EventBinsPos[respondent][clipindex + 1])))
                cmark = TCanvas("c", "c", 1200, 800)
                cmark.cd()
                s = TSpectrum()
                bg = s.Background(hist, 20, "Compton same")
                if type == "phasic":
                    hist.Add(bg,-1)
                elif type == "tonic":
                    hist = bg
                hist.Draw()
                ts = TSpectrum()
                np = ts.Search(hist, sigmapeaksinterval, "noMarkov same nobackground", peakamplitude) / float(len(dataset_index_sub))
                cmark.Update()
                dir = os.path.join(os.getcwd(), './../out/peaks'+name+'/'+EventBinsNames[respondent][clipindex]+'/')
                ensure_dir(dir)
                cmark.SaveAs(dir+'sequencepeaks_'+ respondent + '.png')
                cmark.Close()
                #test = EventBinsNames[respondent][clipindex]
                if('phasic'):
                    peakshist.Fill(EventBinsNames[respondent][clipindex], np)
                else:
                    peakshist.Fill(EventBinsNames[respondent][clipindex], np)
                    #phasichist.Fill(EventBinsNames[respondent][clipindex], np)
                #tonichist.Fill(EventBinsNames[respondent][clipindex], np)
                #hist.Scale(orig_integral/hist.Integral())
        # while (j < nbins):
        #     # create windows of peaks_window size (60 sec)
        #     hist.GetXaxis().SetRange(j, j + peaks_window)
        #     # values found to optimal from above plots
        #     peaskpertime1sig.append(s.Search(hist, 4, "noMarkov nodraw", peakamplitude))
        #     j += peaks_window
        #     indices.append(respondent)
        # values1.extend(peaskpertime1sig)
        hist.Delete()
    return peakshist
    # return pd.Series(values1, index=indices)

def getfirstevent(Events_list,EventBinsNames,respondent):
    for ievt in Events_list:
        for index in range(0,len(EventBinsNames)-1):
            if EventBinsNames[respondent][index] == ievt: return index

#The PCA function is not fully implemented below but work in progress.
def findPrincipalComponent(dataset_index_sub, dataarray):

    #array of TPrincipal objects - with
    fullarray = []
    ## loop over individual "clips"
    for clipindex in range(0, len(EventBinsPos[respondent]) - 1):
       if (clipindex % 2 == 0):
            "find principal component for each sequence?!"
            principal = TPrincipal(len(dataset_index_sub), "D")

            w, h = len(dataarray.dropna(axis=0).loc[dataset_index_sub[0]]),len(dataset_index_sub)
            Matrix = [[0 for x in range(w)] for y in range(h)]
            j=0
            for respondent in dataset_index_sub:
                respdataarray = getrespondent(dataarray, respondent)

                for i in range(0, len(respdataarray) - 1):
                    if respdataarray.index[i] > EventBinsPos[respondent][clipindex] and respdataarray.index[i] < EventBinsPos[respondent][clipindex+1]:
                        Matrix[j][i] = respdataarray.values[i]
                j+=1

            for i in range(0, w):
                row = []
                for j in range(0,h):
                    row.append(Matrix[j][i])
                if Matrix[j][i] != 0:
                    principal.AddRow(array('d',row))
            principal.MakePrincipals()
            print principal.Print("MSEV")
            principal.MakeHistograms("pca"+str(j), "XPDES")
            fullarray.append(principal)

    for respondent in dataset_index_sub:
        respdataarray = getrespondent(dataarray, respondent)

        hist = TH1F(str(respondent) + "hist", str(respondent) + "hist", int((respdataarray.index[-1] * binscale)),
                    respdataarray.index[0], respdataarray.index[-1])

        for i in range(0, len(respdataarray) - 1):
            hist.Fill(respdataarray.index[i], respdataarray.values[i])
        hist.GetXaxis().SetRangeUser(EventBinsPos[respondent][clipindex] + (1 / binscale),
                                     EventBinsPos[respondent][clipindex + 1] - (1 / binscale))

    return fullarray

def meaneda(dataset_index_sub, dataarray,type):
    nbins = (len(dataset_event_names[0]) / 2)
    values0 = TH1F("RMSofEDAperseq", "#sigma_{sequence}", nbins, 0, nbins)
    values1 = TH1F("meanEDAperseq", "#mu_{sequence} #pm #sigma_{sequence}", nbins, 0, nbins)
    values2 = TH1F("meanEDAperseqerr", "#mu_{sequence} #pm sigma_{#mu}", nbins, 0, nbins)
    values1.Sumw2()
    ibin = 0
    labels = sorted(dataset_event_names[0])
    for bin in range(0, len(dataset_event_names[0])):
        if bin % 2 == 0:
            ibin += 1
            values1.GetXaxis().SetBinLabel(ibin, labels[bin])
            values2.GetXaxis().SetBinLabel(ibin, labels[bin])

    for respondent in dataset_index_sub:
        respdataarray = getrespondent(dataarray, respondent)

            # loop over series in timeunit of size: peaks_window
        ## First define peaks per minute for full sequence
        for clipindex in range(0, len(EventBinsPos[respondent]) - 1):
            if (clipindex % 2 == 0):
                hist = TH1F(str(respondent) + str(clipindex) + type +"hist", str(respondent) + str(clipindex) + type + "hist", int(len(respdataarray.values)),respdataarray.values.min(), respdataarray.values.max())
                s = TSpectrum()
                if type == 'phasic':
                    phasichist = TH1F("tmphistphasic", "tmphistphasic",
                                int(len(respdataarray.index)), respdataarray.index.min(), respdataarray.index.max())
                    for i in range(0,len(respdataarray.values)):
                        if respdataarray.index[i] > EventBinsPos[respondent][clipindex]:
                            if respdataarray.index[i] <= EventBinsPos[respondent][clipindex+1]:
                                phasichist.Fill(respdataarray.index[i],respdataarray.values[i])
                    bg = s.Background(phasichist, 20, "Compton same")
                    phasichist.Add(bg, -1)
                    for ibin in range(0,phasichist.GetNbinsX()):
                        hist.Fill(phasichist.GetBinContent(ibin))
                    phasichist.Delete()
                elif type == 'tonic':
                    tonichist = TH1F("tmphisttonic", "tmphisttonic",
                                      int(len(respdataarray.index)), respdataarray.index.min(),
                                      respdataarray.index.max())
                    for i in range(0, len(respdataarray.values)):
                        if respdataarray.index[i] > EventBinsPos[respondent][clipindex]:
                            if respdataarray.index[i] <= EventBinsPos[respondent][clipindex + 1]:
                                tonichist.Fill(respdataarray.index[i], respdataarray.values[i])
                    tonichist = s.Background(tonichist, 20, "Compton same")
                    for ibin in range(0, tonichist.GetNbinsX()):
                        hist.Fill(tonichist.GetBinContent(ibin))
                    tonichist.Delete()
                else:
                    for i in range(0,len(respdataarray.values)):
                        if respdataarray.index[i] > EventBinsPos[respondent][clipindex]:
                            if respdataarray.index[i] <= EventBinsPos[respondent][clipindex+1]:
                                hist.Fill(respdataarray.values[i])

                meanfullrange = hist.GetMean()
                sigmafullrange = hist.GetRMS()
                values0.Fill(EventBinsNames[respondent][clipindex], sigmafullrange)
                values1.Fill(EventBinsNames[respondent][clipindex], meanfullrange)
                values2.Fill(EventBinsNames[respondent][clipindex], meanfullrange)
                #values2.SetBinError(values2.FindBin(EventBinsNames[respondent][clipindex]),values2.GetBinError(values2.FindBin(EventBinsNames[respondent][clipindex]))+ / float(len(dataset_index_sub)))
                hist.Delete()
    for ibin in range(0,values2.GetNbinsX()+1):
        values2.SetBinError(ibin,values0.GetBinContent(ibin))
    return values0,values1,values2

def levene(*args, **kwds):
    """
    Perform Levene test for equal variances.

    The Levene test tests the null hypothesis that all input samples
    are from populations with equal variances.  Levene's test is an
    alternative to Bartlett's test `bartlett` in the case where
    there are significant deviations from normality.

    Parameters
    ----------
    array of samples of measurements:  [[x1,...,xn],[y1,..,yn]]
        The sample data, possibly with different lengths
    center : {'mean', 'median', 'trimmed'}, optional
        Which function of the data to use in the test.  The default
        is 'median'.
    proportiontocut : float, optional
        When `center` is 'trimmed', this gives the proportion of data points
        to cut from each end. (See `scipy.stats.trim_mean`.)
        Default is 0.05.

    Returns
    -------
    statistic : float
        The test statistic.
    pvalue : float
        The p-value for the test.

    Notes
    -----
    Three variations of Levene's test are possible.  The possibilities
    and their recommended usages are:

      * 'median' : Recommended for skewed (non-normal) distributions>
      * 'mean' : Recommended for symmetric, moderate-tailed distributions.
      * 'trimmed' : Recommended for heavy-tailed distributions.

    References
    ----------
    .. [1]  http://www.itl.nist.gov/div898/handbook/eda/section3/eda35a.htm
    .. [2]   Levene, H. (1960). In Contributions to Probability and Statistics:
               Essays in Honor of Harold Hotelling, I. Olkin et al. eds.,
               Stanford University Press, pp. 278-292.
    .. [3]  Brown, M. B. and Forsythe, A. B. (1974), Journal of the American
              Statistical Association, 69, 364-367

    """
    # Handle keyword arguments.
    center = 'median'
    proportiontocut = 0.05
    for kw, value in kwds.items():
        if kw not in ['center', 'proportiontocut']:
            raise TypeError("levene() got an unexpected keyword "
                            "argument '%s'" % kw)
        if kw == 'center':
            center = value
        else:
            proportiontocut = value


    if len(args) > 1:
        raise ValueError("Must supply an iterable list/vector/tuple of samples")
    args=args[0]
    k = len(args)
    if k < 2:
        raise ValueError("Must enter at least two input sample vectors.")
    Ni = zeros(k)
    Yci = zeros(k, 'd')

    if center not in ['mean', 'median', 'trimmed']:
        raise ValueError("Keyword argument <center> must be 'mean', 'median'"
                         + "or 'trimmed'.")

    if center == 'median':
        func = lambda x: np.median(x, axis=0)
    elif center == 'mean':
        func = lambda x: np.mean(x, axis=0)
    else:  # center == 'trimmed'
        args = tuple(scistats.stats.trimboth(np.sort(arg), proportiontocut)
                     for arg in args)
        func = lambda x: np.mean(x, axis=0)

    for j in range(k):
        Ni[j] = len(args[j])
        Yci[j] = func(args[j])
    Ntot = np.sum(Ni, axis=0)

    # compute Zij's
    Zij = [None] * k
    for i in range(k):
        Zij[i] = abs(asarray(args[i]) - Yci[i])

    # compute Zbari
    Zbari = zeros(k, 'd')
    Zbar = 0.0
    for i in range(k):
        Zbari[i] = np.mean(Zij[i], axis=0)
        Zbar += Zbari[i] * Ni[i]

    Zbar /= Ntot
    numer = (Ntot - k) * np.sum(Ni * (Zbari - Zbar)**2, axis=0)

    # compute denom_variance
    dvar = 0.0
    for i in range(k):
        dvar += np.sum((Zij[i] - Zbari[i])**2, axis=0)

    denom = (k - 1.0) * dvar

    W = numer / denom
    pval = scistats.distributions.f.sf(W, k - 1, Ntot - k)  # 1 - cdf
    return W, pval

def f_oneway_custom(*args):
    """
    Performs a 1-way ANOVA.

    The one-way ANOVA tests the null hypothesis that two or more groups have
    the same population mean.  The test is applied to samples from two or
    more groups, possibly with differing sizes.

    Parameters
    ----------
    [sample1, sample2, ...] : list of arrays
        The sample measurements for each group.

    Returns
    -------
    statistic : float
        The computed F-value of the test.
    pvalue : float
        The associated p-value from the F-distribution.

    Notes
    -----
    The ANOVA test has important assumptions that must be satisfied in order
    for the associated p-value to be valid.

    1. The samples are independent.
    2. Each sample is from a normally distributed population.
    3. The population standard deviations of the groups are all equal.  This
       property is known as homoscedasticity.

    If these assumptions are not true for a given set of data, it may still be
    possible to use the Kruskal-Wallis H-test (`scipy.stats.kruskal`) although
    with some loss of power.

    The algorithm is from Heiman[2], pp.394-7.


    References
    ----------
    .. [1] Lowry, Richard.  "Concepts and Applications of Inferential
           Statistics". Chapter 14.
           http://faculty.vassar.edu/lowry/ch14pt1.html

    .. [2] Heiman, G.W.  Research Methods in Statistics. 2002.

    .. [3] McDonald, G. H. "Handbook of Biological Statistics", One-way ANOVA.
           http://http://www.biostathandbook.com/onewayanova.html
    """

    args = [np.asarray(arg, dtype=float) for arg in args[0]]
    # ANOVA on N groups, each in its own array
    num_groups = len(args)
    alldata = np.concatenate(args)
    bign = len(alldata)

    # Determine the mean of the data, and subtract that from all inputs to a
    # variance (via sum_of_sq / sq_of_sum) calculation.  Variance is invariance
    # to a shift in location, and centering all data around zero vastly
    # improves numerical stability.
    offset = alldata.mean()
    alldata -= offset

    sstot = scistats.stats._sum_of_squares(alldata) - (scistats.stats._square_of_sums(alldata) / float(bign))
    ssbn = 0
    for a in args:
        ssbn += scistats.stats._square_of_sums(a - offset) / float(len(a))

    # Naming: variables ending in bn/b are for "between treatments", wn/w are
    # for "within treatments"
    ssbn -= (scistats.stats._square_of_sums(alldata) / float(bign))
    sswn = sstot - ssbn
    dfbn = num_groups - 1
    dfwn = bign - num_groups
    msb = ssbn / float(dfbn)
    msw = sswn / float(dfwn)
    f = msb / msw

    prob = scistats.stats.special.fdtrc(dfbn, dfwn, f)   # equivalent to stats.f.sf

    return f, prob

def ANOVA(dataset_index_sub, dataarray,timewindow):
    """Calculate the intervariance using levenes test due to its robustness over F-test for non-normality in distributions.
    Calculation is performed for each sequence/bin and returned as a two histograms, one with the result per bin and the other with p-value"""
    #glboal_mean_hist = TH1F('glboal_mean_hist', 'glboal_mean_hist', 100, 0, 1)
    # define two histos for the teststatic and p-value results for each sequence
    histtest_median = TH1F("levene_teststatistic_median", "levene_teststatistic_median", int(test_length/timewindow+1), 0, test_length)
    histpval_median = TH1F("levene_pval_median", "levene_pval_median", int(test_length/timewindow+1), 0, test_length)
    histtest_anova = TH1F("anova_teststatistic_anova", "anova_teststatistic_anova", int(test_length/timewindow+1), 0, test_length)
    histpval_anova = TH1F("anova_pval_anova", "anova_pval_anova", int(test_length/timewindow+1), 0, test_length)
    #respondent_bin_dist = TH1F("respondent_bin_dist", "respondent_bin_dist", 20001, 0, 20)

    tasks = []
    binmin = 0
    binsamples = []

    # create multiprocessing pool to allow for faster processing of bins
    #while binmin < int(test_length/timewindow+1):
    #    tasks.append((dataarray, dataset_index_sub, timewindow, binmin,))
    #    binmin += 1
    #pool = multiprocessing.Pool(multiprocessing.cpu_count())

    #print "starting multiprocessing with " + str(multiprocessing.cpu_count()) + " cores"
    # take steps of binsize here and calculate and store in an array to get mean
    for binmin in range(int(test_length/timewindow+1)):
        #binsamples = [pool.apply_async(getSamplesPerTimebin, t) for t in tasks]
        samples = getSamplesPerTimebin(dataarray, dataset_index_sub, timewindow, binmin)

    #print "starting to fill histograms"

    # fill histograms with results of the levene test (F-test alternative) for each bin
    #for binmin in range(int(test_length / timewindow + 1)):
        W, pval = levene(samples, center='median')
        histtest_median.Fill(binmin * timewindow, W)
        histpval_median.Fill(binmin * timewindow, pval)
        #W2, pval2 = levene(samples, center='mean')
        W2, pval2 = f_oneway_custom(samples)
        histtest_anova.Fill(binmin * timewindow, W2)
        histpval_anova.Fill(binmin * timewindow, pval2)
    print "done filling histograms"

    return histtest_median, histpval_median, histtest_anova, histpval_anova


def getSamplesPerTimebin(dataarray, dataset_index_sub, timewindow, binmin):
    samples = []
    for respondent in dataset_index_sub:
        try:
            respdataarray = dataarray.loc[respondent]
        except Exception, e:
            print respondent + "was not in dataarray index - skipping respondent for plots"
            print >> f, respondent + "was not in dataarray index - skipping respondent for plots - error: " + e
            respdataarray = []
        isample = []
        # for clipindex in range(0, len(EventBinsPos[respondent]) - 1):
        #    if (clipindex % 2 == 0):
        start = EventBinsPos[respondent][0] + binmin * timewindow
        end = EventBinsPos[respondent][-1]
        for i in range(0, len(respdataarray.values)):
            if respdataarray.index[i] < start + timewindow * 1000:
                if respdataarray.index[i] >= start:
                    if respdataarray.index[i] < end:
                        isample.append(float(respdataarray.values[i]))
        samples.append(isample)
    return samples

def prep_and_save_hist(thehist, name, title):
    thehist.GetYaxis().SetTitle(title)
    thehist.GetYaxis().SetTitleOffset(1.4)
    thehist.SetStats(0)
    thehist.SetFillColor(4)
    thehist.SetFillStyle(3004)
    thehist.SetLineColor(4)
    thehist.SaveAs('./../out/rootfiles/hist_'+name+'.root')
    # meaneda_full_range.Write()
    c = TCanvas("c", "c", 1200, 800)
    c.cd()
    thehist.Draw("LBAR3")
    c.Update()
    c.SaveAs(os.path.join(os.getcwd(), './../out/results/'+name+'.png'))
    c.SaveAs('./../out/rootfiles/'+name+'.root')

def prep_and_save_hist_plain(thehist, name, title):
    thehist.GetYaxis().SetTitle(title)
    thehist.GetYaxis().SetTitleOffset(1.4)
    thehist.SetStats(False)
    thehist.SetStats(0)
    thehist.SetFillColor(4)
    thehist.SetFillStyle(3004)
    thehist.SetLineColor(4)
    thehist.SaveAs('./../out/rootfiles/hist_'+name+'.root')
    # meaneda_full_range.Write()
    c = TCanvas("c", "c", 1200, 800)
    c.cd()
    thehist.Draw()
    c.Update()
    c.SaveAs(os.path.join(os.getcwd(), './../out/results/' + name + '.png'))
    c.SaveAs(os.path.join(sys.path[0],'../out/rootfiles/'+name+'.root'))