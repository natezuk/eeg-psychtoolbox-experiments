function [corr,resptms] = visual_detect_trial(audio_stim,Fs,visual_stims,targets,wS,varargin)
% This test is meant to deliver an audio stimulus to a participant while the
% participant performs a visual task that consists of button presses in 
% response to a visual target. The audio stimulus is preceded by a click followed by
% a brief period of silence (default 1 s) in order to trigger an Arduino
% to signal the start of the audio presentation in the EEG acquisition 
% hardware (the trigger output from the computer can be somewhat unreliable).
% At the end of this function, the response times are compared to the target 
% times in order to get the number of correct responses.
% Inputs:
% - stim = stimulus audio waveform
% - Fs = sampling rate of the audio (in Hz)
% - audiotrig = stimulus trigger, the value that will be sent to Biosemi to
%       signify which trial it is
% Outputs:
% - corr = number of correct target detections
% - resptms = subject response times
% Bobby Crews & Nate Zuk (2017-2018)

audioprt = []; % audio port
audiotrig = 0; % parallel port trigger for the audio
prllprt = []; % parallel port object, with which to send triggers
sil = 1; % duration of silence before click (in s)
clkamp = 0.6; % magnitude of the click (in V)
clkdur = 1; % duration of the click (in ms)
click_to_stim_time = 1; % amount of time between the click and the stimulus start
scaling_of_texture_square = 10; % scaling of the grid texture relative to the screen
visual_ioi = 0.5; % time between grid presentations
null_texture = {}; % can contain the texture to display on screen before the stimulus starts
detecttol = 0.8; % duration following event at which the keypress is a correct detection (in s)
dbdrop = -35; % amount to drop the stimulus amplitude in dB V to produce a comfortable
    % sound level for stimulus presentation
resptrig = 128; % trigger for subject response

if ~isempty(varargin)
    for n = 2:2:length(varargin)
        eval([varargin{n-1} '=varargin{n};']);
    end
end

% Get the display rectangle on screen for the grid textures
dstRect = [0 0 s1 s2] .* scaling_of_texture_square;
windowRect = Screen('Rect', wS.ptr);
[xCenter, yCenter] = RectCenter(windowRect);
dstRect = CenterRectOnPointd(dstRect, xCenter, yCenter);

if isempty(audioprt),
    error('PsychPortAudio port not specified');
end

%% Present stimulus
% Add click to buffer
clkidx = round(clkdur/1000*Fs); % duration of the click, in indexes
clk = [clkamp*ones(clkidx,1); zeros(Fs*click_to_stim_time-clkidx,1)]; % click followed by 1 second of silence
buf(1) = PsychPortAudio('CreateBuffer',audioprt,ones(2,1)*clk');
% Add stimulus to buffer
audio_stim = adjdb(audio_stim,dbdrop); % set the db V of the stimulus so the sound level is comfortable
buf(2) = PsychPortAudio('CreateBuffer',audioprt,ones(2,1)*audio_stim');

%% Set up audio port schedule with the sound files
% Add stimulus to schedule
PsychPortAudio('UseSchedule',audioprt,2);
PsychPortAudio('AddToSchedule',audioprt,buf(1)); % added click
PsychPortAudio('AddToSchedule',audioprt,buf(2)); % added stimulus

% Send a trigger for the start of the trial (1 second before click)
if ~isempty(prllprt), 
    outputSingleScan(prllprt,dec2binvec(audiotrig,8))
    outputSingleScan(prllprt,dec2binvec(0,8))
end
% After silence, start playing the click and then the stimulus
PsychPortAudio('Start',audioprt,1,GetSecs+sil,1);
%%% Start the visual stimulus here too
strttm = GetSecs+click_to_stim_time; % start time of the stimulus is 1 s after click
keyhold = 0; % flag to check if key is being held down during the trial
resptms = []; % store response times relative to sound start
prtstatus = PsychPortAudio('GetStatus',audioprt);

% Display the null texture, or an empty screen if a null texture isn't
% specified
if isempty(null_texture)
    wS = waitscreen([],wS);
else
    Screen('DrawTextures', wS.ptr, null_texture, [], dstRect, 0, 0); 
    Screen('Flip', wS.ptr);
end
current_stim_display = 0; % no grid is shown yet (before the stimulus starts)

%% During the audio presentation
% While the audio is running (GetStatus:Active is 1)
while prtstatus.Active,
    
    % Display the visual stimuli every visual_ioi
    time_elapsed = GetSecs-strtm;
    % if it is time for a new stimulus
    if ceil(time_elapsed/visual_ioi)~=current_stim_display
        % go to the next visual stimulus
        current_stim_display=ceil(time_elapsed/ioi);
        Screen('DrawTextures', wS.ptr, visual_stims{current_stim_display},...
            [], dstRect, 0, 0);
        Screen('Flip', wS.ptr);
    end
    
    % Check when the subject presses the spacebar
    [keyDown,secs,~] = KbCheck();
    % if the key has been pressed...
    if keyDown,
%             if keyDown==1 && strcmp(KbName(keyCode),'space'),
        % ...check if the key is currently being held down...
        if ~keyhold, %...and if it is...
            % send trigger for keypress
            if ~isempty(prllprt), 
                outputSingleScan(prllprt,dec2binvec(resptrig,8))
                outputSingleScan(prllprt,dec2binvec(0,8))
            end
            resptms = [resptms; secs-strttm-prtstatus.PredictedLatency]; % save the response time
                % adjusted by the start time of the stimulus and
                % the predicted latency of the system
            keyhold=1; % flag that the key is being held down
        end
    else % if the key is not pressed...
        keyhold=0; % ...flag that the key is not being held down
    end
    
    prtstatus = PsychPortAudio('GetStatus',audioprt); % update the state of the audio port
end

%% After audio presentation
% Wait until it ends
PsychPortAudio('Stop',audioprt,1);
% Send a trigger for the end of the trial
if ~isempty(prllprt), 
    outputSingleScan(prllprt,dec2binvec(audiotrig,8))
    outputSingleScan(prllprt,dec2binvec(0,8))
end

%% Compute the number of correct detections
% Compute correct target times
targettms = (targets-1)*visual_ioi; % visual target times, relative to audio start (and visual start)
corr = false(length(targettms),1);
for ii = 1:length(targettms), % for each vibrato time
    % check the number of key presses that were a correct detection of the
    % target time
    detected = sum(resptms>=targettms(ii) & resptms<=targettms(ii)+detecttol);
    % as long as there's one correct keypress, count it as a correct
    % detection
    corr(ii) = detected>0;
end