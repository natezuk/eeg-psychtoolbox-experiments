% Load bdf file and associated results file for SpeechMusicGen_test and
% examine trigger times

addpath('EEGLAB/eeglab14_0_0b/plugins/Biosig3.1.0/biosig/t250_ArtifactPreProcessingQualityControl/');
addpath('EEGLAB/eeglab14_0_0b/plugins/Biosig3.1.0/biosig/t200_FileAccess/');

datapth = 'C:\Users\nzuk\Data\SpVMusRes\testmatlabeeg\';
resnm = 'SpeechMusicGen_res_sbjlprKlv_31-Jan-2018';
bdfnm = 'Testdata3.bdf';

% Load subject responses and results
load([datapth resnm]);

% Load the bdf file
hdr = sopen([datapth bdfnm]);
[fulleeg,hdr] = sread(hdr);
POS = hdr.BDF.Trigger.POS;
TYP = hdr.BDF.Trigger.TYP;

% In the trigger array TYP, the last 7 bits are always set to 1.  Compute
% the value by which to shift the trigger array.
trigadj = binvec2dec([0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1]);

% Plot the triggers
eFs = 512; % sampling rate of bdf file
figure
stem(POS/eFs,TYP-trigadj,'k');
xlabel('Time (s)');
ylabel('Trigger');