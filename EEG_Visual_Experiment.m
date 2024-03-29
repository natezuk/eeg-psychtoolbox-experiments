function results = EEG_Visual_Experiment(varargin)
% Subjects are presented with audio while performing a visual detection task
% Nate Zuk (2019)

% Screen('Preference', 'SkipSyncTests', 1); %This will skip sync tests, should test on pc

% Initial variables
sbj = 'test2';    %subject tag       
stimdir = '/Users/EmilyPrzysinda/Documents/MATLAB/Music_Expectation/Stim/Stim_sel/'; % stimulus directory
visual_ioi = 0.5; % time between grid presentations during the trial
% exmpstimdir = '/Users/EmilyPrzysinda/Documents/MATLAB/Music_Expectation/Classical/All_classical_piano/Stim/ExampleStim/'; % example stimuli directory

Fs = 44100; % should be the same sampling rate as the stimuli (hz)
screens = Screen('Screens'); %Gets the number of screens
scrnnum = max(screens); %once there are multiple screens, select the max value screen (0=home screen)

% Parse varargin (to change initial variables, if desired)
if ~isempty(varargin)
    for n = 2:2:length(varargin)
        eval([varargin{n-1} '=varargin{n};']);
    end
end
    
% Load the list of stimuli
stimlst = dir(stimdir);
% this will be specific to the types of stimuli being used
speechtrknms = stimlst(4:end);

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
objs.start.txt = ['In each trial, you will hear speech or music playing for one minute.\n',...
    'While this is being played, a series of grid patterns will be displayed every 0.5 seconds.\n',...
    'Keep your eyes fixated on the cross in the center and try to blink as little as possible. \n',...
    '\n',...
    'If you''re told to attend the visual stimulus:\n',...
    'Hit the spacebar when you see a grid pattern that is identical to the pattern before\n',...
    'the previous one.\n',...
    '\n',...
    'If you''re told to attend the audio:\n',...
    'Ignore the visual stimulus and pay attention to the speech or music.\n',...
    '\n',...
    'There are ' num2str(length(speechtrknms)) ' trials.\n', ...
    'Press the spacebar when you''re ready.'];
objs.start.active = 1;

% Crosshair shown during trial
objs.cross.type = 'crs';
objs.cross.active = 0;

% Ready text
objs.ready.type = 'dsc';
objs.ready.txt = sprintf('You have completed trial %d/%d\n',...
    'Press any key to start the next trial.',0,0);
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
datafn = ['EEGVisual_res_sbj' sbj '_' date];

%% Setup the trial order
speechtrkorder = randperm(length(speechtrknms)); % randomly rearrange speeech tracks

for jj = 1:length(speechtrkorder) % for each block    
    %% Start the experiment      
    
    % Load the stimulus
    fn=speechtrknms(speechtrkorder(jj)).name;
    [stim Fs]= audioread([stimdir,fn]);
    stim = rampstim(stim,Fs);

    % Show the track names
    disp(speechtrknms{speechtrkorder(jj)});
    
    % Generate the visual stimulus here for the trial
    [grids,targets] = generate_grids(size(stim,1)/Fs,visual_ioi);
    grid_textures = make_grid_screens(wS,grids);
    
    % Create a pre-stimulus display, which is shown before the stimulus
    % starts (and after the click)
    ncols = size(grids,2);
    nrows = size(grids,1);
    null_grid_texture = make_grid_screens(wS,zeros(nrows,ncols));
       
    % run the trial
    [corr,resp] = visual_detect_trial(stim,Fs,grid_textures,targets,...
        wS,null_grid_texture,'audiotrig',speechtrkorder(jj),'visual_ioi',visual_ioi,...
        'audioprt',audioprt,'prllprt',prllprt);
          
    % Save the results
    results.speechtrkorder = speechtrkorder;
    results.speechtrknms = speechtrknms;
    results.corr = corr;
    results.resp = resp;
    save(datafn,'results');
    
    % Display trial number completed and correct responses
    if jj == length(speechtrknms), % if it's the last trial
        objs.end.active = 1;
        [wS,objs] = waitscreen(wS,objs); 
    else
        objs.ready.txt = sprintf(['You have completed trial %d/%d\n',...
            'Press any key to start the next trial.'],jj,length(speechtrkorder));
        objs.ready.active = 1;
        [wS,objs] = waitscreen(wS,objs); 
        objs.ready.active = 0;
    end
end
%% Close PsychToolbox stuff
PsychPortAudio('Close',audioprt);
Screen('Close',wS.ptr);
% Close parallel port
%  clear prllprt ch
end