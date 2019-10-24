function [corr,resptms] = visual_detect_trial(audio_stim,Fs,visual_stim,targettms,stimtrig,varargin)
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
% - targettms = times where then target occurs in the audio
% - stimtrig = stimulus trigger, the value that will be sent to Biosemi to
%       signify which trial it is
% Outputs:
% - corr = number of correct target detections
% - resptms = subject response times
% Bobby Crews & Nate Zuk (2017-2018)
    
wS = []; % window pointer
objs = []; % objects in window
audioprt = []; % audio port
prllprt = []; % parallel port object, with which to send triggers
sil = 1; % duration of silence before stimulus (in s)
clkamp = 0.6; % magnitude of the click (in V)
clkdur = 1; % duration of the click (in ms)
detecttol = 1.5; % duration following event at which the keypress is a correct detection (in s)
dbdrop = -35; % amount to drop the stimulus amplitude in dB V to produce a comfortable
    % sound level for stimulus presentation
resptrig = 128; % trigger for subject response

if ~isempty(varargin)
    for n = 2:2:length(varargin)
        eval([varargin{n-1} '=varargin{n};']);
    end
end

if ~isempty(wS) % (if using psychtoolbox)
    % Reset the screen
    [wS,objs] = gen_screen(objs,wS);
    % Render screen
    Screen('Flip',wS.ptr); 
else
    error('Screen pointer not specified');
end

if ~isempty(audioprt),
    %% Present stimulus
    % Add click to buffer
    clkidx = round(clkdur/1000*Fs); % duration of the click, in indexes
    clk = [clkamp*ones(clkidx,1); zeros(Fs-clkidx,1)]; % click followed by 1 second of silence
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
        outputSingleScan(prllprt,dec2binvec(stimtrig,8))
        outputSingleScan(prllprt,dec2binvec(0,8))
    end
    % After silence, start playing the click and then the stimulus
    PsychPortAudio('Start',audioprt,1,GetSecs+sil,1);
    %%% Start the visual stimulus here too
    strttm = GetSecs+1; % start time of the stimulus is 1 s after click
    keyhold = 0; % flag to check if key is being held down during the trial
    resptms = []; % store response times relative to sound start
    prtstatus = PsychPortAudio('GetStatus',audioprt);
    
    %% During the audio presentation
    % While the audio is running (GetStatus:Active is 1)
    while prtstatus.Active,
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
    % Send a trigger for the end of the trial (1 second before click)
    if ~isempty(prllprt), 
        outputSingleScan(prllprt,dec2binvec(0,8))
        outputSingleScan(prllprt,dec2binvec(0,8))
    end
else
    error('PsychPortAudio port not specified');
end

%% Compute the number of correct detections
corr = false(length(targettms),1);
for ii = 1:length(targettms), % for each vibrato time
    % check the number of key presses that were a correct detection of the
    % target time
    detected = sum(resptms>=targettms(ii) & resptms<=targettms(ii)+detecttol);
    % as long as there's one correct keypress, count it as a correct
    % detection
    corr(ii) = detected>0;
end