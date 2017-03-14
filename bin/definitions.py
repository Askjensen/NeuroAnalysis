# coding=utf-8

# output to file:
f = open('../out/output.txt', 'w')
# f = open(os.path.join(os.getcwd(), '../out/output.txt'), 'w')
# ASSUMPTIONS #
# timesteps of 32ms so 48*32/1000 s = 1,536 s
# 5 * 1000 / 32 = 156,25 ~ 156
# 10 * 1000 / 32 ~ 312
test_length = 4000 #the total length in seconds of the test including baselinevideos 5000s/60s/min =
rolling_average_window = 156  # number of data points to include within a rolling avg window (p=0.3 pearson's corr between original EDA data and normalized phasic data)
peaks_window = 60  # 60 bins of 1 sec. = 1 min windows for peak finding
timewindow = 10 #sek bins for peak finding and mean EDA
binscale = 1./1000. #32. / 60000.  # 1 sek/bin
sigmapeaksinterval = 3 # the significance required above background for peaks
peakamplitude = 0.01 #relative amplitude to highest peak required
endbuffer = 0

#### datafile tag definitions - these match naming scheme of Biometric Software Suite output-files
sync_pos = 'position'
eda_data = 'EDA'
pupil_data = ['PupilLeft','PupilRight']
event_data = 'tag__info_StudioEventData'
delimiter = ';'

## add constant to ensure no negative EDA values - if not distributions normalised to unity can flip EDA signal sign for some reagions.
addconstant = True

# List of events in the data - these are the names defined in tobii and should be given as a list
Events_list = ['f490aca9-f27c-4328-a2a4-f1ab259c481c.avi','baseline.avi','Bankerot_1.avi','Bankerot_2.avi']
#Events_list = ['Baseline.avi']

# A list of sequences to be directly compared

Comparison_list = {'1':['Bankerot_1.avi','Bankerot_2.avi']}
#                    '2':['03_Absalon_og_Valdemar.avi','03_Absalon_og_Valdemar.avi'],
#                    '3':['04_Absalon_og_Valdemar.avi','04_Absalon_og_Valdemar.avi'],
#                    '4':['05_Absalon_og_Valdemar.avi','05_Absalon_og_Valdemar.avi'],
#                    '5':['06_Saxo.avi','06_Saxo.avi'],
#                    '6':['07_Ingeborg.avi','07_Ingeborg.avi'],
# #                    '7':['08_Ingeborg.avi','08_Ingeborg.avi'],
# #                    '7': ['09_Ingeborg.avi', '09_Ingeborg.avi'],
# #                    '7': ['10_Ingeborg.avi', '10_Ingeborg.avi'],
# #                    '7': ['11_Valdemar_Sejr.avi', '11_Valdemar_Sejr.avi'],
# #                    '7': ['12_Valdemar_Sejr.avi', '12_Valdemar_Sejr.avi'],
# #                    '7': ['13_Valdemar_Sejr.avi', '13_Valdemar_Sejr.avi'],
# #                    '7': ['14_Ærkebispen_af_Lund.avi', '14_Ærkebispen_af_Lund.avi'],
# #                    '7': ['15_Erik_Klipping.avi', '15_Erik_Klipping.avi'],

#Path to the directory containing datafiles (output ASCII txt-files from Biometric Software Suite)
folder_path = 'C:/Data/Bankerot/lette/'
#folder_path = 'C:/Data/HistorienOmDK/alder/over 50/'
#folder_path = 'C:/Data/HistorienOmDK/alder/under 40/'
#folder_path = 'C:/Data/HistorienOmDK/women/'
#folder_path = 'C:/Data/HistorienOmDK/males/'
#folder_path = 'C:/Data/HistorienOmDK/alleinklfejl/'
#folder_path = 'C:/Data/HistorienOmDK/rest/'
#folder_path = 'C:/Data/HistorienOmDK/test/'

filename_ext = '.txt'
index1_names = list()
index2_syncpos = list()
eda_values_list = list()
pupil_diameter_values_list = list()
filelist = list()
dataset_index = list()
dataset_event_names = list()
dataset_event_datapoints = list()
EventBinsPos = {}
EventBinsNames = {}