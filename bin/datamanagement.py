import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from bin.Functions import dataextract
from bin.definitions import f, endbuffer, sync_pos, eda_data, pupil_data, event_data, delimiter, Events_list, \
    folder_path, filename_ext, dataset_index, dataset_event_names, EventBinsPos, dataset_event_datapoints, \
    EventBinsNames, index1_names, index2_syncpos, eda_values_list, pupil_diameter_values_list, filelist


def event_hz_pointers(files, events, event_data_sub, sync_pos_sub):
    fhandle = open(folder_path + files, 'r')

    n = 0
    event_temp_list = list()
    event_hz_markers = list()
    previous = ""

    for p in fhandle:

        p1 = p.split(delimiter)
        value_to_add = 0
        if n == 0:
            label_lookup_event = p1.index(event_data_sub)
            label_lookup_position = p1.index(sync_pos_sub)

        if n != 0:
            # if len(p1) < label_lookup_event:
            #     print >> f, "Error. Unknown error occured when slicing events from dataset"
            #     print >> f, "A discrepancy between lenght of dataset and event position. Event position exceeded the lenght of dataset."
            #     print >> f, "Lenght of dataset: ", len(p1)
            #     print >> f, "Event label position: ", label_lookup_event
            #     print >> f, "track back error halts script from continuing."
            if len(p1) >= label_lookup_event:
                #print len(p1)
                #print label_lookup_event
                #print label_lookup_position
                #test = 0
                if p1[label_lookup_event] in events:
                    event_temp_list.append(p1[label_lookup_event])

                    try:
                        value_to_add = int((p1[label_lookup_position]))
                    except:
                        print >> f, "Error. Blank cell or string data detected where integer data was expected."
                        print >> f, "This will result in failure to slice events. Dataset should be corrected before continuing"
                        print >> f, "Stop script execution and remove or correct dataset."
                    #Add a buffer to the end of the sequence to catch EDA reactions with delay.
                    if(p1[label_lookup_event]==previous):
                        event_hz_markers.append(value_to_add+endbuffer)
                    else:
                        event_hz_markers.append(value_to_add)
                    previous = p1[label_lookup_event]
        n = n + 1

    fhandle.close()

    event_temp_list1, event_hz_markers = (list(x) for x in
                                         zip(*sorted(zip(event_temp_list, event_hz_markers), key=lambda pair: pair[0])))

    return (event_temp_list, event_hz_markers)


def eventmeans(name, sub_event_names, sub_event_datapoints, df):
    temp_mean = list()
    temp_std = list()

    number_of_events = len(sub_event_names)
    n = 0

    while n < number_of_events:
        get_slice_start = int(sub_event_datapoints[n])
        get_slice_stop = int(sub_event_datapoints[n + 1])

        event_arousal_mean = df.loc[name].loc[get_slice_start:get_slice_stop].mean()
        event_arousal_std = df.loc[name].loc[get_slice_start:get_slice_stop].std()

        temp_mean.append(event_arousal_mean.values)
        temp_std.append(event_arousal_std.values)

        n = n + 2

    return (temp_mean, temp_std)


def suplabel(axis, label, label_prop=None,
             labelpad=5,
             ha='center', va='center'):
    ''' Add super ylabel or xlabel to the figure
    Similar to matplotlib.suptitle
    axis       - string: "x" or "y"
    label      - string
    label_prop - keyword dictionary for Text
    labelpad   - padding from the axis (default: 5)
    ha         - horizontal alignment (default: "center")
    va         - vertical alignment (default: "center")
    '''
    fig = plt.gcf()
    xmin = []
    ymin = []
    for ax in fig.axes:
        xmin.append(ax.get_position().xmin)
        ymin.append(ax.get_position().ymin)
    xmin, ymin = min(xmin), min(ymin)
    dpi = fig.dpi
    if axis.lower() == "y":
        rotation = 90.
        x = xmin - float(labelpad) / dpi
        y = 0.5
    elif axis.lower() == 'x':
        rotation = 0.
        x = 0.5
        y = ymin - float(labelpad) / dpi
    else:
        raise Exception("Unexpected axis: x or y")
    if label_prop is None:
        label_prop = dict()
    plt.text(x, y, label, rotation=rotation,
             transform=fig.transFigure,
             ha=ha, va=va,
             **label_prop)


def read_data():
    # Walk through folder creating file list
    for roots, dirs, files in os.walk(folder_path):
        for files_to_compute in files:
            if files_to_compute.endswith(filename_ext):
                filelist.append(files_to_compute)
    print >> f, "Number of datasets to compute: ", len(filelist)
    # Extract data from files
    for masterfile in filelist:
        print >> f, "Processing dataset: ", masterfile
        index_list = dataextract(masterfile, sync_pos)
        eda_data_list = dataextract(masterfile, eda_data)
        pupil_data_array = []
        pupil_data_list = []
        for i in range(0, len(pupil_data)):
            pupil_data_array.append(dataextract(masterfile, pupil_data[i]))
        for i in range(0, len(pupil_data_array[0])):
            pupil_data_list.append((pupil_data_array[0][i] + pupil_data_array[1][i]) / 2.0)

        print masterfile
        event_names, event_datapoints = event_hz_pointers(masterfile, Events_list, event_data,
                                                          sync_pos)  # syntax: file, list of events to look for, data tag in data set describing events, sync positions, delimiter
        id_name = (masterfile.split('.'))[0]

        if len(index_list) != len(eda_data_list):
            "WARNING! in file", id_name, " a difference between index lenght and value lenght has been encounted."
        if len(index_list) != len(pupil_data_list):
            "WARNING! in file", id_name, " a difference between index lenght and pupil diameter lenght has been encounted."
        n1 = 0
        id_name_list = list()

        while (len(index_list)) > n1:
            id_name_list.append(id_name)
            n1 = n1 + 1

        # Clean data focusing solely on defined events. Other data = NaN.

        xn = 0
        xn1 = 0
        xn2 = 0

        eda_data1 = list()
        pupil_data1 = list()
        event_hz_markers1 = event_datapoints
        event_hz_markers1.sort()

        for xp in index_list:

            if len(event_hz_markers1) != xn and event_hz_markers1[xn] == xp:
                # Switching - when in event data array add data otherwise add numpy.nan
                xn = xn + 1
                xn1 = 0 if xn1 else 1

            if xn1 == 0:
                eda_data1.append(np.nan)  # Fill NaN when outside event boundaries
                pupil_data1.append(np.nan)

            elif xn1 == 1:
                eda_data1.append(eda_data_list[xn2])
                pupil_data1.append(pupil_data_list[xn2])

            xn2 = xn2 + 1

        # All data are extended to three lists, essentially creating a MultiIndex (index1+2) and referring values
        index1_names.extend(id_name_list)
        index2_syncpos.extend(index_list)
        eda_values_list.extend(eda_data1)
        pupil_diameter_values_list.extend(pupil_data1)

        # This list keeps track of all data sets in order to compute event means
        dataset_index.append(id_name)
        dataset_event_names.append(event_names)
        dataset_event_datapoints.append(event_datapoints)
        EventBinsPos[id_name] = event_datapoints
        EventBinsNames[id_name] = event_names


def create_dataframes():
    #global eda_data_series, pupil_data_series
    print >> f, "Creating dataframe..."
    # Make MultiIndex a Tuple zipping two lists together
    index3 = list(zip(index1_names, index2_syncpos))
    # Create 2D MultiIndex
    index4 = pd.MultiIndex.from_tuples(index3, names=['Names', 'Syncpos'])
    # Parse the Index to create a new dataframe named DF
    df = pd.DataFrame(index=index4)
    # Create a new series with the original data
    eda_data_series = pd.Series(eda_values_list, index=index4)
    # resultdf = eda_data_series.to_frame()
    # print resultdf.head
    # print resultdf.dropna(axis=0).head
    pupil_data_series = pd.Series(pupil_diameter_values_list, index=index4)
    return eda_data_series.to_frame().dropna(axis=0),pupil_data_series.to_frame().dropna(axis=0)