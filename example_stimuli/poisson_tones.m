% Create a click train with Poisson distribution timing
% (based on Maddox & Lee, 2018)
% Present tones randomly sampled between 800-1250 Hz

addpath('..');

click_rate = 2; % average click rate in s
stim_dur = 30; % duration of the click train
ntr = 20; % number of trials
Fs = 48000;
tone_dur = 50; % duration in ms
tone_freq = [800 1250]; % tone frequencies, in Hz
rtime = 5; % duration of onset/offset ramps for tones

t_tone = (0:round(tone_dur/1000*Fs)-1)/Fs; % time array for the tone

stim_path = 'tones/';
ct_path = 'tone_times/';

for n = 1:ntr
    % Generate the randomized timings of the clicks
    ioi = exprnd(1/click_rate,ceil(click_rate*stim_dur*2),1);
    % remove any iois less than the tone duration
    ioi(ioi<tone_dur/1000) = [];
    onset = cumsum(ioi);
    onset(onset>stim_dur-tone_dur/1000) = []; % remove any click times longer than the stimulus duration
    
    % Initialize the stimulus array
    stim = zeros(stim_dur*Fs,1);
    for ii = 1:length(onset) % for each onset

        % Create the tone
        f = rand*diff(tone_freq)+tone_freq(1); % randomly select a tone frequency
        tone = sin(2*pi*t_tone*f);
        % apply a 5 ms onset/offset ramp
        tone = rampstim(tone,Fs,'rtime',rtime);
        
        % Place the clicks at the correct tone at the correct time 
        idx = round(onset(ii)*Fs)+(0:length(tone)-1);
        stim(idx) = tone;

    end

    % Save the audio
    fl = sprintf('tones_%02g',n);
    audiowrite([stim_path fl '.wav'],stim,Fs);
    disp(fl);

    % Save the click times
    writematrix(onset,[stim_path ct_path fl '.csv']);
end