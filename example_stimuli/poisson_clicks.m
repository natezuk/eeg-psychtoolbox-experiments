% Create a click train with Poisson distribution timing
% (based on Maddox & Lee, 2018)

click_rate = 2; % average click rate in s
stim_dur = 30; % duration of the click train
ntr = 20; % number of trials
Fs = 48000;
stim_path = 'clicks/';
ct_path = 'click_times/';

for n = 1:ntr
    % Generate the randomized timings of the clicks
    ioi = exprnd(1/click_rate,ceil(click_rate*stim_dur*2),1);
    click_times = cumsum(ioi);
    click_times(click_times>stim_dur) = []; % remove any click times longer than the stimulus duration
    
    % Create the click (100 us duration rarefaction click)
    click_dur = 1e-4;
    click_idx = 0:ceil(click_dur*Fs)-1;
    
    % Place the clicks at the correct times 
    stim = zeros(stim_dur*Fs,1);
    for ii = 1:length(click_times)
        idx = round(click_times(ii)*Fs);
        stim(click_idx+idx) = -1; % rarefaction click
    end

    % Save the audio
    fl = sprintf('clicks_%02g',n);
    audiowrite([stim_path fl '.wav'],stim,Fs);
    disp(fl);

    % Save the click times
    writematrix(click_times,[stim_path ct_path fl '.csv']);
end