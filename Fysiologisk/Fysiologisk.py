#!/usr/bin/env python
# coding=utf-8
######## Ask Emil Løvschall-Jensen, July 2016 ## DR Audience Research ##############
# The script reads in datafiles combining GSR (scin-conductance) data and eye-tracking (pupil-dillation) and calculates agregated distributions.
# Functions are defined first in the script and the 'main' afterwards. The CERN statistical frameworks ROOT and spectral analysis algoritms
# from TSpectrum are used to model data. TSpectrum is used to find the phasic and Tonic component of data as well as phasic peak positions.
# Normalised distributions of these are combined. By default, the background is removed before deconvolution in TSpectrum.
# We specify the option "nobackground" to not remove the background as we instead model it and manually substract it.

#Input is defined in the bin/defintion.py file

# The scripts give the following output under the out folder:
#folder:        content:
#peaksEDA       Individual respondent distrubtions of Scin-conductance data for all sequnces and for full overview
#peaksPD        Individual respondent distrubtions of Pupil-dillation data for all sequnces and for full overview
#phasicpeaks    The same plots as in the above folders but from the phasic component of data alone, so peaks are found from phasic data.
#results        numberofpeaks* plots for number of peaks per sequence and tonic and phasic component
#results        timedistributionsof* plots for number of peaks VS time for each sequence in comparison_list and tonic and phasic component
#output.txt     log-file with all outputs from running the script - including som standards parameters for distributions - not very tidy!

## Define which plots to generate - set to True for plots to be made
#Create distributions of threshold impact on peaks - use only on baseline data with known peaks.
import cProfile

testthtresholds = False
#Create distributions for each respondent with EDA data
peakseda = False #done
#Create distributions for each respondent with pupil-dillation data
peakspd = False #done
#Create distributions for each respondent with normalized phasic component of data
phasic = False #done

### Create full range plots of eda, phasic, tonic and peak values in each timebin with size
## defined in the definition file (timewindow, default: 10 seconds).
rawedapeaks = False
rawpdpeaks = False

#plot test-statistic and p-value for ANOVA and levene test of variations
doAnova = False

#plot mean values for sequences:
meanraw = False #done
meanphasic = True
meantonic = False
#should the mean be calculated per sequence or for each timebin. True for Sequences, False for timebins
meanFromSequences = False
#do further overvieplots: needs revision, not working and disabled
dooverview = False

import bin.definitions
from bin.Functions import testsettings, createplotsFullRange, phasic_component, tonic_component, normalize_series, meaninterval, \
    npeaks, meaneda, prep_and_save_hist, prep_and_save_hist_plain, findPrincipalComponent, ANOVA
from bin.datamanagement import suplabel, read_data, create_dataframes
from bin.definitions import f, rolling_average_window, timewindow, Comparison_list
import matplotlib.pyplot as plt
# from bin.BH import *
# ROOT is a statistical tool from CERN - it can be downloaded here: https://root.cern.ch/content/release-53436
# release: root_v5.34.36.win32.vc12.exe (windows release, but available for MAC and most linux distributions)) (requires VC++ 12.0)
# import ROOT

# sys.path.append('/Applications/root_v5.34.36/lib')
from ROOT import TFile, TCanvas, gStyle, TPrincipal, TSpectrum

#plot styling can be defined here:
#import myrootstyle
# -*- coding: utf-8 -*-


#
from matplotlib.font_manager import FontProperties
import os

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

#profile the functioncalls and timeusage of the program? This is quite time-consuming. Only enabled functions will be tested.
doProfile=False


def runmain():
    ### PROGRAM ###
    print >> f, "Script to perform automated data analysis of GSR and eye-tracking data."
    print >> f, "Ask E. Loevschall-Jensen 2016"
    print >> f, "---------------------------------------------------------------------------"
    print >> f, "Data from the specified dir will now be read"


    # Initialize global lists

    read_data()

    eda_data_series,pupil_data_series = create_dataframes()


    print >> f, eda_data_series
    # Tonic/phasic computations
    # Window size = number of datapoint!

    print >> f, "Computing number of peaks in time-window"

    rootfile = TFile(os.path.join(os.getcwd(), '../out/rootfiles/histos.root'), "RECREATE")

    if (testthtresholds):
        testsettings(bin.definitions.dataset_index, eda_data_series)
    #PCAs = findPrincipalComponent(bin.definitions.dataset_index,eda_data_series)
    #test = 0

    if doAnova:
        anova_hist_median,anova_pval_hist_median,anova_hist_mean,anova_pval_hist_mean = ANOVA(bin.definitions.dataset_index, eda_data_series,timewindow)
        prep_and_save_hist(anova_hist_median,'levene_teststatistic_median','levene teststatistic using median for ANOVA per timebin')
        prep_and_save_hist(anova_pval_hist_median,'levene_pval_median','levene p-value using median for ANOVA per timebin')
        prep_and_save_hist(anova_hist_mean,'anova_teststatistic','ANOVA teststatistic per timebin')
        prep_and_save_hist(anova_pval_hist_mean,'anova_pval','ANOVA p-value per timebin')

    if(peakseda):
        ################ find peaks, phasic and tonic distributions in the data ######################
        ################ all distributions are normalised to  integral 1 #############################
        print 'peaksperminute_full_range'
        peaksperminute_full_range = npeaks(bin.definitions.dataset_index, eda_data_series, 'EDA','raw')
        prep_and_save_hist(peaksperminute_full_range,'EDApeaksperseq_full_range','Peaks in EDA per sequence')
        # #peaksperminute_full_range.GetXaxis().SetTitle('Time [minutes]')
        # peaksperminute_full_range.GetYaxis().SetTitle('Peaks per sequence')
        # peaksperminute_full_range.GetYaxis().SetTitleOffset(1.4)
        # peaksperminute_full_range.SetStats(False)
        # peaksperminute_full_range.SetFillStyle(3002)
        # peaksperminute_full_range.Write()
        # c = TCanvas("c", "c", 1200, 800)
        # c.cd()
        # peaksperminute_full_range.Draw()
        # c.Update()
        # c.SaveAs(os.path.join(os.getcwd(), './../out/results/EDApeaksperseq_full_range.png'))

    if(peakspd):
        print 'peaksperminute_full_range_PD'
        peaksperminute_full_range_PD = npeaks(bin.definitions.dataset_index, pupil_data_series, 'PD','raw')
        prep_and_save_hist(peaksperminute_full_range_PD,'PDpeaksperseq_full_range','Peaks in PD per sequence')
        #
        # #peaksperminute_full_range_PD.GetXaxis().SetTitle('Time [minutes]')
        # peaksperminute_full_range_PD.GetYaxis().SetTitle('Peaks per sequence')
        # peaksperminute_full_range_PD.GetYaxis().SetTitleOffset(1.4)
        # peaksperminute_full_range_PD.SetStats(False)
        # peaksperminute_full_range_PD.SetFillStyle(3002)
        # peaksperminute_full_range_PD.Write()
        # c = TCanvas("c", "c", 1200, 800)
        # c.cd()
        # peaksperminute_full_range_PD.Draw()
        # c.Update()
        # c.SaveAs(os.path.join(os.getcwd(), './../out/results/PDpeaksperseq_full_range.png'))



    fontP = FontProperties()
    fontP.set_size('small')


    ### Create comparision plots for eda and/or PD data for sequences in Comparison_list
    if(rawedapeaks):
        # plot specific events in same hist
        #npeaksseqtotal= npeaksspecific(dataset_index, eda_data_series, Comparison_list,'rawEDApeaks')
        #npeaksseqminuts = npeaksspecificminutes(dataset_index, eda_data_series, Comparison_list,timewindow,'rawEDApeaks')
        print 'createplotsFullRange(raweda)'
        npeaksfullrange = createplotsFullRange(bin.definitions.dataset_index, eda_data_series, timewindow, 'rawEDApeaks')

    if(rawpdpeaks):
        #npeaksseqtotal= npeaksspecific(dataset_index, pupil_data_series, Comparison_list,'rawPDpeaks')
        #npeaksseqminuts = npeaksspecificminutes(dataset_index, pupil_data_series, Comparison_list,timewindow,'rawPDpeaks')
        print 'createplotsFullRange(rawpd)'
        npeaksfullrange = createplotsFullRange(bin.definitions.dataset_index, pupil_data_series, timewindow, 'rawPDpeaks')


    ###########################################################################
    ###########################################################################

    # peaksinwindow_series1 = npeaks(dataset_index, eda_data_series)



    # a plot peaks in the full datarange in 1 minute bins summed over all respondents.
    # peaksperminute_full_range = TH1F("peaksperminute_full_range", "Peaks per minute with 4 #sigma significance",
    #                                 len(peaksinwindow_series1), 0, len(peaksinwindow_series1))

    # plot peaks for the full datarange and for events defined in Comparison_list
    # for i in range(0, len(peaksinwindow_series1)):
    #    peaksperminute_full_range.Fill(i, peaksinwindow_series1.values[i])
    # eda_data_series.loc['02418'].index[-1] = 2.629.585
    # for event in Comparison_list:
    #    for pos in range(int(EventBinsDict['start_'+event]/(1000. * 60.)),int(EventBinsDict['end_'+event]/(1000. * 60.))):
    #        int(EventBinsDict['start_MF vaccine med case.avi'] / (1000. * 60.))

    # For specific events of interest in same histogram
    # peaksperminute_compared = TH1F("peaksperminute_compared", "Peaks per minute with 4 #sigma significance",
    #                               len(peaksinwindow_series1), 0, len(peaksinwindow_series1))


    ####################### find phasic component of EDA ######################
    ###########################################################################
    # timescale for plots
    timescaling = 1 / (1000. * 1)  #

    if dooverview and False:
        print >> f, "Computing phasic component of data..." ############this is the old definiton of time-avaraging!!!!
        phasic_data_series = phasic_component(bin.definitions.dataset_index, eda_data_series, rolling_average_window)
        tonic_data_series = tonic_component(bin.definitions.dataset_index, eda_data_series, rolling_average_window)
        # print >> f, phasic_data_series
        # t=phasic_data_series.index
        for respondent in bin.definitions.dataset_index:
            # plot raw EDA for first case:
            tmpraw = eda_data_series.loc[respondent]
            tmptonic = tonic_data_series.loc[respondent]
            # Two subplots, the axes array is 1-d
            sp1 = plt.subplot(2, 1, 1)
            plt.plot(tmpraw.index * timescaling, tmpraw.values, label='Raw EDA')
            plt.plot(tmptonic.index * timescaling, tmptonic.values, label='tonic (time avaraged EDA)')
            sp1.set_title(respondent)
            # plt.ylabel('Scin conductance (EDA)')
            plt.legend(bbox_to_anchor=(1.1, 1.2), prop=fontP)
            # phasic part
            plt.subplot(2, 1, 2)
            tmpphasic = phasic_data_series.loc[respondent]
            plt.plot(tmpphasic.index * timescaling, tmpphasic.values,
                     label='phasic rest after subtraction of tonic time avarage')
            plt.legend(bbox_to_anchor=(1.1, 0.15), prop=fontP)
            plt.grid(True)
            plt.xlabel('time [s]')
            suplabel('y', 'Skin conductance (EDA)')
            plt.savefig("./../out/respondents/" + respondent + ".png")
            plt.close()
            # plt.show()

    # Normalize series output
    print >> f, "Normalizing phasic data..."
    raw_normalized = normalize_series(bin.definitions.dataset_index, eda_data_series)

    if (meanraw):
        ######################### mean EDA on normalized EDA #########################
        ##############################################################################
        #find mean eda for each sequence
        print 'raw: sigmaofeda_full_range,meaneda_full_range,meaneda_full_range_err'
        sigmaofeda_full_range,meaneda_full_range,meaneda_full_range_err = meaneda(bin.definitions.dataset_index, eda_data_series,'raw',meanFromSequences) #raw_normalized
        name = "sequence" if meanFromSequences else str(timewindow)+" sec. bin"
        prep_and_save_hist(meaneda_full_range, 'meaneda_full_range', '#mu_{'+name+'} #pm  1 #sigma_{'+name+'}')
        prep_and_save_hist(meaneda_full_range_err, 'meaneda_full_range_err', '#mu_{'+name+'} #pm 1 #sigma_{#mu}')
        prep_and_save_hist_plain(sigmaofeda_full_range, 'sigmaofeda_full_range', '#sigma_{#mu}')


        # find mean per interval for each sequence - histogram as for peaks
        meanedainterval = meaninterval(bin.definitions.dataset_index, raw_normalized, Comparison_list, timewindow, 'rawmean')
        #meanedainterval.GetYaxis().SetTitle('Mean EDA per timeinterval')
        #meanedainterval.GetYaxis().SetTitleOffset(1.4)
        #meanedainterval.SetStats(False)
        #meaneda_full_range.Write()
        #c = TCanvas("c", "c", 1200, 800)
        #c.cd()
        #meanedainterval.Draw()
        #c.Update()
        #c.SaveAs(os.path.join(os.getcwd(), './../out/results/meaneda_intervals.png'))
        ###########################################################################
        ###########################################################################

    if (meantonic):
        ######################### mean EDA on normalized EDA #########################
        ##############################################################################
        #find mean eda for each sequence
        name = "sequence" if meanFromSequences else str(timewindow)+" sec. bin"
        print 'tonic: sigmaofeda_full_range,meaneda_full_range,meaneda_full_range_err'
        sigmaofeda_full_range,meaneda_full_range,meaneda_full_range_err = meaneda(bin.definitions.dataset_index, eda_data_series,'tonic',meanFromSequences)
        prep_and_save_hist(meaneda_full_range, 'meanedatonic_full_range', '#mu_{'+name+'} #pm  1 #sigma_{'+name+'}')
        prep_and_save_hist(meaneda_full_range_err, 'meanedatonic_full_range_err', '#mu_{'+name+'} #pm  1 #sigma_{#mu}')
        prep_and_save_hist_plain(sigmaofeda_full_range, 'sigmaoftoniceda_full_range', '#sigma_{#mu} tonic')

    if(meanphasic): ############this is the old definiton of time-avaraging!!!!
        #phasic_data_series = phasic_component(bin.definitions.dataset_index, eda_data_series, rolling_average_window)
        #phasic_normalized = normalize_series(bin.definitions.dataset_index, phasic_data_series)
        ################ mean EDA based on phasic component of data ###############
        ###########################################################################
        #find mean eda for each sequence of phasic data
        name = "sequence" if meanFromSequences else str(timewindow)+" sec. bin"
        print 'phasic: sigmaofeda_full_range,meaneda_full_range,meaneda_full_range_err'
        sigmaofeda_full_range,meanphasic_full_range,meanphasic_full_range_err = meaneda(bin.definitions.dataset_index, eda_data_series,'phasic',meanFromSequences)
        prep_and_save_hist(meanphasic_full_range, 'meanphasic_full_range', 'Mean Phasic per '+name)
        prep_and_save_hist(meanphasic_full_range_err, 'meanphasic_full_range_err', 'Mean Phasic per '+name)
        prep_and_save_hist_plain(sigmaofeda_full_range, 'sigmaofphasiceda_full_range', 'Standard deviation of Phasic part of EDA per '+name)

        # find phasic mean per interval for each sequence - histogram as for peaks
        #todo: make this with true phasic if interval comparisons are needed:
        # meanphasicinterval = meaninterval(bin.definitions.dataset_index, phasic_normalized, Comparison_list, timewindow, 'phasicmean')

        #eanphasicinterval.GetYaxis().SetTitle('Mean Phasic per interval')
        #meanphasicinterval.GetYaxis().SetTitleOffset(1.4)
        #meanphasicinterval.SetStats(False)
        #meaneda_full_range.Write()
        #c = TCanvas("c", "c", 1200, 800)
        #c.cd()
        #meanphasicinterval.Draw()
        #c.Update()
        #c.SaveAs(os.path.join(os.getcwd(), './../out/results/meanphasic_intervals.png'))
        ###########################################################################
        ###########################################################################


    if(phasic and peakseda): ############this is the old definiton of time-avaraging!!!!
        ############################### peaks phasic ##############################
        ###########################################################################
        ##find and plot histograms for peaks in phasic data
        print 'npeaks oldphasic'
        npeaksseqphasic= npeaks(bin.definitions.dataset_index, eda_data_series, 'phasicEDA','phasic')
        prep_and_save_hist(npeaksseqphasic,'npeaksseqphasic','Peaks per sequence')

        ###########################################################################
        ###########################################################################



    if dooverview and False:
        phasic_data_series = phasic_component(bin.definitions.dataset_index, eda_data_series, rolling_average_window)
        phasic_normalized = normalize_series(bin.definitions.dataset_index, phasic_data_series)
        tonic_normalized = normalize_series(bin.definitions.dataset_index, tonic_data_series)
        ######################### test of individuals #########################
        # plot normalized for short range
        for respondent in bin.definitions.dataset_index:
            # plot raw EDA for first case:
            tmpraw = raw_normalized.dropna(axis=0).loc[respondent]
            tmptonic = tonic_normalized.dropna(axis=0).loc[respondent]
            # Two subplots, the axes array is 1-d
            sp1 = plt.subplot(2, 1, 1)
            plt.plot(tmpraw.index[1000:3000] * timescaling, tmpraw.values[1000:3000], label='Normalized EDA')
            plt.plot((tmptonic.index[1000:3000]) * timescaling, tmptonic.values[1000:3000],
                     label='Normalized tonic (time avaraged EDA)')
            sp1.set_title(respondent)
            # plt.ylabel('Scin conductance (EDA)')
            plt.legend(bbox_to_anchor=(1.1, 1.2), prop=fontP)
            # phasic part
            plt.subplot(2, 1, 2)
            tmpphasic = phasic_normalized.dropna(axis=0).loc[respondent]
            plt.plot(tmpphasic.index[1000:3000] * timescaling, tmpphasic.values[1000:3000], label='Normalized phasic EDA')
            plt.legend(bbox_to_anchor=(1.1, 0.15),   prop=fontP)
            plt.grid(True)
            plt.xlabel('time [s]')
            suplabel('y', 'Normalized skin conductance (EDA)')
            plt.savefig("./../out/respondents/" + respondent + "_normalized.png")
            plt.close()
            # plt.show()

    #
    # if(True):
    #     # add column to existing DataFrame:
    #     df[df.index[1], 'EDA'] = phasic_normalized
    #
    #     # Data validation: Test for equal number of events per dataset,
    #     print >> f, "Validating events for each dataset... Number of events and their names must be exactly the same across all datasets"
    #
    #     n = 0
    #     data_events_validating = len(dataset_index[0])
    #     while (len(dataset_index) - 1) > n:
    #         if (len(dataset_event_names[n])) != (len(dataset_event_names[n + 1])):
    #             print >> f, "DATA FAILURE... Unequal number of events detected!! Means may not be valid!"
    #
    #         n = n + 1
    #
    #     # Initalize new dataframe - dataframe for means #
    #     df1 = pd.DataFrame(index=dataset_event_names[0][0:(len(dataset_event_names[0])):2])
    #
    #     n = 0
    #
    #     for p in dataset_index:
    #         event_means, event_std = eventmeans(p, dataset_event_names[n], dataset_event_datapoints[n], df)
    #         n = n + 1
    #
    #         df1[p] = event_means
    #
    #     # Output event with arousal mean and standard deviation
    #     print >> f, "Output: event name, arousal mean for event and standard deviation for mean."
    #
    #     # asktodo - plot this:
    #     for p in dataset_event_names[0][0:(len(dataset_event_names[0])):2]:
    #         print >> f, p + "," + str((df1.loc[p].mean())[0]) + "," + str(df1.loc[p].std())
    #
    #     # Repeated Measures ANOVA
    #     print >> f, "Performing a repeated measures ANOVA..."
    #
    #     df_anova = pt.DataFrame()
    #
    #     headers = namedtuple('headers', ['subject', 'event', 'mean'])
    #
    #     an = 0
    #
    #     for p in dataset_event_names[0][0:(len(dataset_event_names[0])):2]:
    #
    #         for p1 in dataset_index:
    #             p2 = float(df1.loc[p].loc[p1])
    #             df_anova.insert(headers(p1, an, p2)._asdict())
    #
    #         an = an + 1
    #
    #     anova = df_anova.anova(dv='mean', sub='subject', wfactors=['event'])
    #
    #     print >> f, (anova)
    #
    #     # Concluding the script doing a T-test to test for sig. diff. means between two predefined means
    #     for ikey in Comparison_list.keys():
    #         t_test_list1 = list()
    #         t_test_list2 = list()
    #         print >> f, "Testing for significance between: " + Comparison_list[ikey][0] + " and " + Comparison_list[ikey][1]
    #         t_temp_list1 = df1.loc[ Comparison_list[ikey][0]].values
    #         t_temp_list2 = df1.loc[ Comparison_list[ikey][1]].values
    #
    #         for p in t_temp_list1:
    #             t_test_list1.append(float(p))
    #
    #         for p1 in t_temp_list2:
    #             t_test_list2.append(float(p1))
    #
    #         print >> f, ttest_ind(t_test_list1, t_test_list2)

    rootfile.Write()
    rootfile.Close()

###main###

if doProfile:
    cProfile.run('runmain()', 'profile.result','cumtime')
else:
    runmain()


