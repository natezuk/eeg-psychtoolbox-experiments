% Test vibrato_detect_trial, make sure that it can detect spacebar presses
% as the audio from PsychPortAudio is playing
% Nate Zuk (2018)

% Make a stimulus consisting of regularly timed beeps
addpath('teststimuli');
dur = 20; % duration of stimulus (in s)
Fs = 44100; % sampling frequency (in Hz)
% 120 BPM, 440 Hz, 100 ms pips, randomly select phase
[stim,piptms] = tonepips(20,440,0.1,rand(1)-0.5,dur,Fs);

%% Initialize audio port
InitializePsychSound();
% Check available sound devices
dev = PsychPortAudio('GetDevices');
% Open audio port
prt = PsychPortAudio('Open',[],[],0,Fs,1);

%% Create the screen
objs.instr.type = 'dsc';
objs.instr.txt = 'Press the spacebar each time you hear a beep';
objs.instr.active = 1;
[wS,objs] = gen_screen(objs,[],'dsp',0);

%% Start the trial
[corr,resptms] = vibrato_detect_trial(stim,Fs,piptms,'wS',wS,'objs',objs,'prt',prt);
% Display number of correct detections
fprintf('%d/%d correctly detected\n',sum(corr),length(piptms));

%% Close Psych stuff
PsychPortAudio('Close',prt);
Screen('Close',wS.ptr);