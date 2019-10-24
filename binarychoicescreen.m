function [resp,wS,objs] = binarychoicescreen(wS,objs,keys)
% When using PsychToolbox, update the screen and wait for the user to press
% one of two keys to specify a binary choice. The active text should state 
% what the choice is. 'keys' is a 2-element cell array containing the
% strings for the keys corresponding to 0 and 1 responses respectively
[wS,objs] = gen_screen(objs,wS);
Screen('Flip',wS.ptr);
while 1,
    [keyDown,~,keyCode] = KbCheck();
    % if the key has been pressed...
    if keyDown,
        % ...check which key it is...
        if strcmp(KbName(keyCode),keys{1}), %...and if it is...
            resp = false; % first key is false
            KbReleaseWait(); % wait until the key is released
            return
        elseif strcmp(KbName(keyCode),keys{2}),
            resp = true; % second key is true
            KbReleaseWait(); % wait unti lthe key is released
            return
        end
    end
end