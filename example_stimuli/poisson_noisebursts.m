% Create a click train with Poisson distribution timing
% (based on Maddox & Lee, 2018)
% Present tones randomly sampled between 800-1250 Hz

addpath('..');

click_rate = 2; % average click rate in Hz
stim_dur = 30; % duration of the click train
ntr = 20; % number of trials
Fs = 48000;
noise_dur = 50; % duration in ms
rtime = 5; % duration of onset/offset ramps for tones
amp = 0.1; % RMS amplitude of the noise burst

% t_noise = (0:round(noise_dur/1000*Fs)-1)/Fs; % time array for the tone

stim_path = 'noise/';
ct_path = 'noise_times/';

for n = 1:ntr
    % Generate the randomized timings of the clicks
    ioi = exprnd(1/click_rate,ceil(click_rate*stim_dur*2),1);
    % remove any iois less than the tone duration
    ioi(ioi<noise_dur/1000) = [];
    onset = cumsum(ioi);
    onset(onset>stim_dur-noise_dur/1000) = []; % remove any click times longer than the stimulus duration
    
    % Initialize the stimulus array
    stim = zeros(stim_dur*Fs,1);
    for ii = 1:length(onset) % for each onset

        % Create the noise burst
        noise = randn(round(noise_dur/1000*Fs),1)*amp;
        % apply a 5 ms onset/offset ramp
        noise = rampstim(noise,Fs,'rtime',rtime);
        
        % Place the clicks at the correct tone at the correct time 
        idx = round(onset(ii)*Fs)+(0:length(noise)-1);
        stim(idx) = noise;

    end

    % Save the audio
    fl = sprintf('noise_%02g',n);
    audiowrite([stim_path fl '.wav'],stim,Fs);
    disp(fl);

    % Save the click times
    writematrix(onset,[stim_path ct_path fl '.csv']);
end