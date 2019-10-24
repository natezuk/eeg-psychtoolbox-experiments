function testkeypress()
% Display the key name when the keyboard is pressed

while 1,
    [keydown,~,keycode] = KbCheck();
    if keydown,
%         disp(KbName(keycode));
        if strcmp(KbName(keycode),'right'),
            disp('true'); 
        else
            disp('false');
        end
        break
    end
end