function [resp,wS,objs] = ratingscalechoice(wS,objs)
% When using PsychToolbox, update the screen and wait for the user to chose 1-5 on a rating scale
%The active text should state 
% what the choice is. 'keys' is a 2-element cell array containing the
% strings for the keys corresponding to 0 and 1 responses respectively
[wS,objs] = gen_screen(objs,wS);
Screen('Flip',wS.ptr);
 while 1,
        [keyDown,~,keyCode,~] = KbCheck();
        if keyDown,
            if strcmp(KbName(keyCode),'1!'),
                resp = 1;
                KbReleaseWait(); % wait until the key is released
                return
            elseif strcmp(KbName(keyCode),'2@'),
                resp = 2;
                KbReleaseWait(); % wait until the key is released
                return
            elseif strcmp(KbName(keyCode),'3#'),
                resp = 3;
                KbReleaseWait(); % wait until the key is released
                return
            elseif strcmp(KbName(keyCode),'4$'),
                resp = 4;
                KbReleaseWait(); % wait until the key is released
                return
            elseif strcmp(KbName(keyCode),'5%'),
                resp = 5;
                KbReleaseWait(); % wait until the key is released
                return
            elseif strcmp(KbName(keyCode),'esc'),
                resp = NaN;
                KbReleaseWait(); % wait until the key is released
                return
            end
        end

    end