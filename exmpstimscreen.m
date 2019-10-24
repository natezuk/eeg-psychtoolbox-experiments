function [wS,objs] = exmpstimscreen(wS,objs,prt,keys,stims,varargin)
% Using psychtoolbox, create the start screen displaying instruction text
% and allowing subjects to play example sounds with key presses. Inputs are 'obj',
% the structure containing on-screen objects, and 'prt' for the
% PsychPortAudio audio port.
% Nate Zuk (2018)

% Initial variables
exmppth = '';
dbdrop = -35; % amount to drop the stimulus amplitude in dB V to produce a comfortable
    % sound level for stimulus presentation

if ~isempty(varargin),
    for n=2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

% Show the new screen
[wS,objs] = gen_screen(objs,wS);
Screen('Flip',wS.ptr);

% Get the buttons on screen
objnms = fieldnames(objs);
btnnms = {};
for n = 1:length(objnms),
    eval(['typ = objs.' objnms{n} '.type;']);
    if strcmp(typ,'btn'), % if it's a button
        btnnms = [btnnms objnms(n)]; % save it in the list of buttons
    end
end

% Wait for user responses
while 1,
    % Wait for the user to press the left or right arrow key to play an
    % example stimulus
%     [clx,cly,btn] = GetMouse(wS.ptr);
    [keyDown,~,keyCode] = KbCheck();
    if keyDown,
        % Check if the subject pressed a key corresponding to an example
        % stimulus
        keychk = strcmp(KbName(keyCode),keys);
        if sum(keychk)==1,
            % load the corresponding sound file...
            [snd,~] = audioread([stims{keychk}]);
            snd = adjdb(snd,dbdrop); % set db of stimulus
            buf = PsychPortAudio('CreateBuffer',prt,ones(2,1)*snd');
            % ...and play it
            PsychPortAudio('UseSchedule',prt,2);
            [success,~] = PsychPortAudio('AddToSchedule',prt,buf);
            if ~success, 
                warning('Could not add stimulus to audio buffer.');
            else,
                PsychPortAudio('Start',prt,1);
                PsychPortAudio('Stop',prt,1);
            end
            KbReleaseWait(); % wait for the key to be released
        end
        % Check if the subject pressed spacebar
        if strcmp(KbName(keyCode),'space'),
            return
        end
    end
end