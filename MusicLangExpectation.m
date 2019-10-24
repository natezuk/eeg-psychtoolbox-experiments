function results = MusicLangExpectation(varargin)
% Subjects are presented with a stimulus of speech or music and are asked
% to decide if the stimulus was scrambled

Screen('Preference', 'SkipSyncTests', 1); %This will skip sync tests, should test on pc

% Initial variables
sbj = 'test2';    %subjct tag       
stimdir = '/Users/EmilyPrzysinda/Documents/MATLAB/Music_Expectation/Stim/Stim_sel/'; % stimulus directory
exmpstimdir = '/Users/EmilyPrzysinda/Documents/MATLAB/Music_Expectation/Classical/All_classical_piano/Stim/ExampleStim/'; % example stimuli directory


Fs = 44100; % should be the same sampling rate as the stimuli
screens = Screen('Screens'); %Gets the number of screens (will only be 1 for now)
scrnnum = max(screens); %once there are multiple screens, select the max value screen (0=home screen)
exmpkeys = {'right','down','left'}; % key responses to play different example stimuli
bnrykeys = {'right','left'}; % key responses to make a binary choice
exmpnms = {'musicwobbleexmp.wav','noiseexmp.wav','speechwobbleexmp.wav'};

% Parse varargin (to change initial variables, if desired)
if ~isempty(varargin)
    for n = 2:2:length(varargin)
        eval([varargin{n-1} '=varargin{n};']);
    end
end


%Loading the audio files
    
% Load the list of stimuli
stimlst = dir(stimdir);
speechtrknms = stimlst(4:end);
%order = randperm(length(speechtrknms));

%% Initialize audio port
InitializePsychSound();
% % Check available sound devices
% dev = PsychPortAudio('GetDevices');
% Open audio port
%audioprt = PsychPortAudio('Open',6,1,1,Fs,2);
% Using Sound Blaster Audio speaker output in sound booth (NZ, 1/29/2018)

audioprt = PsychPortAudio('Open',[],[],1,Fs,2); 

%% Initialize the parallel port interface
prllprt = [];
% ch = addDigitalChannel(prllprt,'Dev1', 'Port2/Line0:7', 'OutputOnly'); % specify the channels to use
%     % (port 2 is connected with the USB receiver)
% outputSingleScan(prllprt, [0 0 0 0 0 0 0 0]) % reset parallel port to 0

%% Setup the screen
% Start text
objs.start.type = 'dsc';
objs.start.txt = ['In each trial, you will hear one song that is 2-3 minutes in length.\n',...
    'While you are listening to the song, keep your eyes fixated on the cross in the center. \n',...
    'Try to blink as little as possible while the trial is playing.\n',...
    '\n',...
    'In between each song you will be asked to rate your familiarity with the song \n',...
    'and how much you liked it on a scale from 1-5\n',...
    'There are ' num2str(length(speechtrknms)) ' trials.\n', ...
    'Press the spacebar when you''re ready.'];
objs.start.active = 1;

% Crosshair shown during trial
objs.cross.type = 'crs';
objs.cross.active = 0;

% Recognition text
objs.recog.type = 'dsc';
objs.recog.txt = ['Rate how familiar you were with the song on a scale from 1 to 5:\n',...
    '\n',...
    'Use the number keys on the keyboard for the following rating scale\n',...
    '1   2   3   4   5\n\n',...
    '1 - I have never heard this song before in my life\n',...
    '2 - I may have heard this song or a song like it\n',...
    '3 - I have heard this song at least once before\n',...
    '4 - I have heard this song multiple times before and have a general sense of how it sounds\n',...
    '5 - I have heard this song many times and am very familiar with it  2  3  4  5'];

objs.recog.active = 0;

% Like text
objs.like.type = 'dsc';
objs.like.txt = ['Rate how you liked listening to this song\n',...
    'on a scale from 1 being extremely dislike to 5 being extremely like',...
    '\n',...
    'Use the numbers on the keyboard for the rating scale:',...
    '1  2  3  4  5'];
objs.like.active = 0;

% Ready text
objs.ready.type = 'dsc';
objs.ready.txt = ['Press the space bar when you are ready for the\n',...
    '\n', 'next trial'];
objs.ready.active = 0;

% End text
objs.end.type = 'dsc';
objs.end.txt = ['Congrats, you have finished the experiment!\n',...
    '\n',...
    'Press any key to exit.'];
objs.end.active = 0;

%% Make the screen and display initial instructions
[wS,objs] = gen_screen(objs,[],'dsp',scrnnum);
[wS,objs] = exmpstimscreen(wS,objs,audioprt,exmpkeys,exmpnms,'exmppth',exmpstimdir);
objs.start.active = 0; % turn off instruction text

% Reset random number generator
rng('shuffle')

%% Create filename
if isempty(sbj),
    sval = randi([65,90],1,6); % randomly pick 6 uppercase letters (ASCII values)
    sbj = char(sval); % convert the values to ASCII
end
datafn = ['SpeechMusicGen_res_sbj' sbj '_' date];

%% Setup the trial order
speechtrkorder = randperm(length(speechtrknms)); % randomly rearrange speeech tracks
% trkorder = d.results.trkorder(14:end);
% wblnum = randi(3,length(trkorder),1); % randomly select 1-3 for each trial to determine number of wobbles

recogresp = zeros(length(speechtrknms),1); % logical array to store familiarity responses
likeresp = zeros(length(speechtrknms),1); % logical array to store like responses

for jj = 1:length(speechtrkorder) % for each block    
    %% Start the experiment      
    
    % Load the stimulus
    fn=speechtrknms(speechtrkorder(jj)).name;
    [stim Fs]= audioread([stimdir,fn]);
    stim = rampstim(stim,Fs);

    % Show the track names
    %disp(speechtrknms{speechtrkorder(jj)});
       
    % run the trial
    objs.cross.active = 1; % show the crosshair
    vibrato_detect_trial(stim,Fs,[],...
        'wS',wS,'objs',objs,'audioprt',audioprt,'prllprt',prllprt,'stimtrig',speechtrknms(speechtrkorder(jj)).name);
    objs.cross.active = 0; % turn off the crosshair
              
    % Recognize the song?
    objs.recog.active = 1;
    [recogresp(jj,1),wS,objs] = ratingscalechoice(wS,objs);
    objs.recog.active = 0;
    
    % Like the song?
    objs.like.active = 1;
    [likeresp(jj,1),wS,objs] = ratingscalechoice(wS,objs);
    objs.like.active = 0;
    
    % Ready? 
    objs.ready.active = 1;
    [wS,objs] = waitscreen(wS,objs);
    objs.ready.active = 0;
          
    % Save the results
    results.recogresp = recogresp;
    results.likeresp = likeresp;
    results.speechtrkorder = speechtrkorder;
    results.speechtrknms = speechtrknms;
    save(datafn,'results');
    
    % Display trial number completed and correct responses
    if jj == length(speechtrknms), % if it's the last trial
        objs.break.txt = sprintf('Congrats, you have finished the experiment!\n',...
            '\n',...
            'Press any key to exit.');
        objs.end.active = 1;
        [wS,objs] = waitscreen(wS,objs); 
    end
end
%% Close PsychToolbox stuff
PsychPortAudio('Close',audioprt);
Screen('Close',wS.ptr);
% Close parallel port
%  clear prllprt ch
end