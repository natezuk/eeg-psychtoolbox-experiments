function nshv_v8b(subj_id, acq_no, cb_order, which_task_to_do, varargin)
% function nshv_v8b(subj_id, acq_no, cb_order, which_task_to_do, varargin)
%
% <which_task_to_do>: 'auditory' || 'visual'
% 
% started from: nshv_v8a_demoSpatialNoSound_bitFlip_proportionTime.m
%
% modified from v8a:
%   -- record the visual responses
%   -- record the auditory responses
%   -- don't display the visual task during the NULL conditions
% 
% similar to nshv_v5, but playing the mixture of the FGs with BG-shaped noise
%     rather than the BGs by themselves
% 
% high variation natural sounds experiments (mixing with different
%   levels of background noise)
% 
% <subj_id> is a string that you have the corresponding stim for
%     -- IT MUST BE OF FORMAT: 'sub%d' -- string where first three characters are 'sub'
%     -- the last index is used to index into the stim_cb_orders
% <acq_no> doesn't do anything, but is used to label files
% <cb_order> determines the counterbalance order
% <varargin> can be anything -- then you're in debugging mode
%
% CREATED: 06.09.2017
% [ modified from nshv_v5,
%   which was modifed from nshv_v3b,
%     which was modified from nshv_v3, 
%       which was modified from nshv_v2, 
%         which was modified from my fall's speech pilot code, 
%           which was from my 2013 textures code, 
%             which sure looks like it was modified from my mgh code...]
% 
% LAST MODIFIED: 08.30.2018
%
% CONDITIONS:
%   54 different wavs:
%       27 unique "sources" x 2 conditions 
%           [ just FG, FG+BG]
%   -- but then two different run types: either visual attention or auditory attention
%
% there are 18 possible cb orders for each subject (and they are different across subjects)
%
% runs for: 347.490056 seconds
%   TR: 3.270 seconds
%   ips = 106 [including the TRs we're tossing away at beginning 
%               -- sparse acq. protocol doesn't have disc. acqs.]
%   5:47.49 min:sec

% check if debugging
if ismember('debug',varargin) || ismember('DEBUG',varargin) % do debug
  PsychDebugWindowConfiguration(0,0.4);
  dbstop if error
  ShowCursor();
else
  % PsychDebugWindowConfiguration(0,1.0);
  dbclear all
  clear Screen
  HideCursor();
end

if ismember('demo',varargin) || ismember('DEMO',varargin) % demo so save accordingly
  is_demo = true;
else
  is_demo = false;
end

which_task_to_do = lower(which_task_to_do);
if ~ismember(which_task_to_do,{'auditory','visual'})
  error('<which_task_to_do> has to be ''visual'' or ''auditory'' (note that it''s not case sensitive\n  you entered: %s\n',...
              which_task_to_do);
end

% set up for run
prefs = setPreferences(subj_id, acq_no, cb_order);
[win, prefs] = genericInitialize(prefs, is_demo, which_task_to_do);
prefs = initializeAudio(prefs);


% get stim info 
[event_info,prefs] = getEventInfo(prefs, win);
event_info = preloadAudio(event_info,win,prefs);


% normalize onsets for run to the trigger time
run_onset_absolute = waitForTrigger(win, prefs, which_task_to_do);
event_info.absolute_onsets = event_info.relative_onsets + run_onset_absolute;


% info string for the top of the print out screen... 
info_str = 'onset\tdur\ttruedur\tcond#\tcond\tstimName\tkey\trt';
fprintf(['\n\n\n' info_str '\n\n']);


% start the timer...
tic

for event_idx = 1:length(event_info.conds)
  
  actual_absolute_onset = GetSecs();
  event_info = presentEvent(event_info,event_idx,win,prefs,which_task_to_do);
  actual_absolute_offset = GetSecs();
  actual_event_length = actual_absolute_offset - actual_absolute_onset;
  
  % output info about event
  cond_str = prefs.conditions{event_info.conds(event_idx)+1};
  if event_info.conds(event_idx)==0 % fixation
    event_str = sprintf('%.2f\t%.3f\t%.3f\t%d\t%s\t%s\t\n', ...
                         event_info.relative_onsets(event_idx), event_info.durs(event_idx), actual_event_length, ...
                         event_info.conds(event_idx), cond_str,event_info.clip_names{event_idx});
  else
    switch which_task_to_do
      case 'auditory'
        key_press = event_info.key_pressed_auditory(event_idx,:);
        rt = event_info.rt_auditory(event_idx);
      case 'visual'
        key_press = event_info.key_pressed_visual(event_idx,:);
        rt = event_info.rt_visual(event_idx);
    end
    event_str = sprintf('%.2f\t%.3f\t%.3f\t%d\t%s\t%s\t\t%s\t%.3f\n', ...
                         event_info.relative_onsets(event_idx), event_info.durs(event_idx), actual_event_length, ...
                         event_info.conds(event_idx), cond_str,event_info.clip_names{event_idx}, ...
                         num2str(key_press), rt);
  end

  fprintf('%s',event_str);
  fprintf(prefs.dataFid,event_str);
end
toc

doEndOfRunHouseKeeping(prefs,event_info);

end


function p = setPreferences(subj_id, acq_no, cb_order)

if ~strcmp(subj_id(1:3),'sub')
  error('wrong format for the subj_id -- must start with ''sub'' and then have a number');
end

p = struct;
p.subj_id = subj_id;
p.sub_no = str2double( subj_id(4:end) ); % used to access into run orders
p.acq_no = acq_no;
p.cb_order = cb_order;
p.monitor = max(Screen('Screens'));
p.bg_color = ones(1,3)*80; % a nice grey
p.fg_color = [1 1 1]*255;

p.incorrect_color = [1 0 0] * 255;
p.correct_color = [0 1 0] * 255;

p.trigger = [KbName('+') KbName('=+')];
p.response_keys = [ KbName('1!') KbName('1')]; 
% p.responseKeys = [ KbName('1!') KbName('2@') KbName('3#') KbName('4$')];
% p.responseKeys = [KbName('1') KbName('2') KbName('3') KbName('4')];
p.quit_keys = [KbName('q')]; %#ok<NBRAK>
p.size_dot = 20;
p.text_size = 48;
exp_dir = '/Users/mcdermottscanner/Desktop/alex/nshv-v8b-attention/';
p.data_dir = [exp_dir filesep 'data' filesep];
p.stim_dir = [exp_dir filesep 'stim/final-stim/' subj_id filesep];
p.stim_cb_order_fpath = [exp_dir filesep 'counterbalance-orders/stim_cb_orders.mat'];

p.time = struct();

% TR: 2.96s; 2 cyc: 1.48s; 3 cyc: 0.986s; 4 cyc: 0.74s
% p.time.tr = 2.96;

% TR: 3.27s; 2 cyc: 1.635s, 3 cyc: 1.09s; 4 cyc: 0.8175s
p.time.tr = 3.27;

% p.time.n_cycles_per_tr = 2; 
% p.time.n_cycles_per_tr = 3;  
p.time.n_cycles_per_tr = 4; % chose to do four for the v8b experiment

% p.time.prop_vis_stim = 0.7; p.time.prop_vis_isi = 0.3;
p.time.prop_vis_stim = 0.85; p.time.prop_vis_isi = 0.15;

p.time.vis_trial_length = (p.time.tr / p.time.n_cycles_per_tr); 
p.time.vis_stim_length = p.time.vis_trial_length * p.time.prop_vis_stim;
p.time.vis_isi = p.time.vis_trial_length * p.time.prop_vis_isi;

% obviously the three proportions need to sum to 1
assert(1 == p.time.prop_vis_stim + p.time.prop_vis_isi);
assert(abs(p.time.tr-p.time.n_cycles_per_tr * p.time.vis_trial_length)<10^-10); % within a tolerance

% p.n_back = 2; % 1-back is just detect repeats
p.n_back = 1; % 1-back is just detect repeats
p.probability_of_repeat = 0.3; % as a result, mean = 0.2283 repeated, std = 0.0226 across 100 samples
% p.probability_of_repeat = 0.5;
% p.probability_of_repeat = 1.0; % for debuggin'

% p.n_bits_flipped = 1;
p.n_bits_flipped = 2;

p.audio_task_prop_visual_randomly_pseudocorrect = 0.9; % and randomly incorrect 1 - p

p.design = struct();
p.design.n_trials_per_block = 3;
p.time.acq_time = 0.870; % previous: 0.560

p.design.n_trs_dummy_at_beginning = 4;
p.design.n_trs_null_beg = 3;
p.design.n_trs_null_end = 6;

% params for debuggin':
% p.design.n_trs_dummy_at_beginning = 0;
% p.design.n_trs_null_beg = 0;
% p.design.n_trs_null_end = 6;

% GRID PARAMETERS
p.grid = struct();

p.grid.n_filled_in = 6;
p.grid.scaling_of_texture_square = 10; 
p.grid.off_square_color = p.bg_color;
p.grid.on_square_color = [204 0 204]; % magenta

% want a 4x4 grid
p.grid.grid_width_in_squares = 4; 
p.grid.grid_height_in_squares = 4;

% parameterize the square itself
p.grid.square_edge = 9; 
p.grid.grid_width = 1;


p.freq = 44100; % sampling rate of stimuli
p.stim_ext = 'wav';

p.conditions = {'NULL'};
for ii = 0:(54-1) % zero indexing FTW; NOTE: **this is 54 hard-coded!!
  p.conditions{end+1} = sprintf('stim%02.f',ii);
end

p.audio_device = 2; % check which devices are there with PsychPortAudio('GetDevices')
end


function [win,p] = genericInitialize(p,is_demo,which_task_to_do)

commandwindow();
% ListenChar(2); % don't let characters output to commandwindow
% HideCursor();
addpath('/Users/mcdermottscanner/Desktop/alex/sam-lib/');
% rng('shuffle'); % seed the random number generator with the time

Screen('Preference', 'VisualDebugLevel', 1); % make the initial screen black instead of white
Screen('Preference','SkipSyncTests',1); % AK added 2018.07.26 because failing tests: WARNING: Couldn't compute a reliable estimate of monitor refresh interval! Trouble with VBL syncing?!?
win = Screen('OpenWindow',p.monitor,p.bg_color);
Screen('FillRect',win,p.bg_color);
Screen('Flip',win);
Screen('TextSize',win,p.text_size);
% version_string = PsychtoolboxVersion();
% version_number = regexp(version_string,'[0-9]*\.[0-9]*\.[0-9]*','match');
% drawCenteredText(win,sprintf('Initializing Psychtoolbox...\n\nVersion number: %s',version_number{1}),p.fg_color);
Screen('Flip',win);

tail_str = '';
if is_demo; tail_str = [tail_str '_DEMO']; end
tail_str = [tail_str '_' upper(which_task_to_do) '.txt'];

p.data_fpath = [p.data_dir filesep 'nshv-v8b_' p.subj_id '_acqNo' num2str(p.acq_no) '_cb' num2str(p.cb_order) tail_str]; 

if exist(p.data_fpath,'file')
  p.data_fpath = [p.data_fpath(1:end-4) '.' datestr(now,30) '.txt'];
end
p.dataFid = fopen(p.data_fpath,'w');

end


function p = initializeAudio(p)

InitializePsychSound();

mode = 1; % 1 is just playback
reqlatencyclass = 1;  % how aggressively to pursue deterministic latency;
                      % 1 is default; can increase up to 4, which is greedy 
                      % and aggressive
% freq = 20000; % requested playback latency in Hz

% PsychPortAudio('Open' [, deviceid][, mode = 1 for playback][, reqlatencyclass][, freq][, channels][, buffersize][, suggestedLatency][, selectchannels][, specialFlags=0]);
p.pa_handle = PsychPortAudio('Open',p.audio_device-1,mode,reqlatencyclass,p.freq); % note ppa uses zero-idxing
% if need lower-latency configurations check out help InitializePsychSound();
end


function [e,p] = getEventInfo(p,w)
% conditions: the 90 wav files 
% NOTE: cond==0 null -- no sound and no grid

% each row is a counter balance order
% not including intial or final fixation (which is the length of a block)
% generated these from sam's scripts -- removed ones with fixation at end or beginning

tmp = load(p.stim_cb_order_fpath);
cb_orders = tmp.stim_cb_orders;
clear tmp

% **change this if first stim is delayed...
initial_fix_length = p.time.acq_time + p.time.tr * p.design.n_trs_dummy_at_beginning + p.time.tr * p.design.n_trs_null_beg; 
final_fix_length = p.time.tr * p.design.n_trs_null_end; % go through the last TR

% note that adding one to counteract the 0-idxing; also: (n_sub, n_runs_per_sub, n_blocks_per_run)
cb_order = squeeze(cb_orders(p.sub_no+1, p.cb_order, :)); 
n_blocks = length(cb_order);

% get all clip names
stims_st = struct;
d = dir(p.stim_dir);
all_clip_names = {d.name}; 
for cond_idx = 2:length(p.conditions) % first condition is null
  condition = p.conditions{cond_idx};
  % appending this underscore/dot is vital to prevent me from being an idiot. as is the error below.
  is_clip = cellfun(@(c)(~isempty(strfind(c,[condition '.']))), all_clip_names); %#ok<STREMP>
  clip_name = all_clip_names(is_clip);
  
  if length(clip_name)>1
    error('<clip_name>: glob is matching more than one file. problems!');
  end
  stims_st.(condition) = clip_name{1}; % fucking stupid matlab cell arrays
end

% generate event_info and assign clip names

% initialize event_info and do first fixation 
e = struct;
e.relative_onsets = 0;
e.conds = 0; % because NULL
e.durs = initial_fix_length;
e.clip_names = {''};
e.is_quieter_clip = nan;

% make one of the clips quieter: 10^(-7/20)

for block_idx = 1:n_blocks
  condition_no = cb_order(block_idx);
  
  % refactoring this to assign stim blocks and ISIs to each of the null block components
  if condition_no==0
    quieter_trial_idx = nan; % nan == x --> returns false always
  else
    quieter_trial_idx = randi(p.design.n_trials_per_block-1) + 1; % don't have the first be quieter...
  end

  for stim_idx = 1:p.design.n_trials_per_block
    % stim
    e.relative_onsets(end+1) = e.relative_onsets(end) + e.durs(end);
    e.conds(end+1) = condition_no;
    e.durs(end+1) = p.time.tr;

    if condition_no==0
      e.clip_names{end+1} = '';
    else
      e.clip_names{end+1} = stims_st.(p.conditions{condition_no+1}); 
    end
    
    % determine whether clip is quieter or not
    % NOTE: for condition_no==0 --> no clip is quieter because there is no clip
    if stim_idx==quieter_trial_idx
        e.is_quieter_clip(end+1) = true;
    else
        e.is_quieter_clip(end+1) = false;
    end
  end
end

% end of run fixation, too
e.relative_onsets(end+1) = e.relative_onsets(end) + e.durs(end); 
e.conds(end+1) = 0;
e.durs(end+1) = final_fix_length;
e.clip_names{end+1} = '';
e.is_quieter_clip(end+1) = nan;

% create the key-pressed field
n_events = length(e.conds);
e.key_pressed_visual = nan * ones(n_events,p.time.n_cycles_per_tr);
e.rt_visual = nan * ones(n_events,p.time.n_cycles_per_tr); 

e.key_pressed_auditory = nan * ones(n_events,1);
e.rt_auditory = nan * ones(n_events,1);


% actually generate the grids (textures) that you want
%   in the order you want with the appropriate repeats
recent_grids = nan * ones(p.grid.grid_height_in_squares, ...
                          p.grid.grid_width_in_squares, ...
                          p.n_back);
e.grids = cell(n_events,p.time.n_cycles_per_tr);
e.grid_textures = cell(n_events,p.time.n_cycles_per_tr);
e.grid_is_repeated = nan * ones(n_events,p.time.n_cycles_per_tr);

% generate initial idxs randomly
n_els_in_grid = p.grid.grid_width_in_squares*p.grid.grid_height_in_squares;
% randomly sample without replacement
previous_grid_idxs = randsample(n_els_in_grid,p.grid.n_filled_in,false); 

is_previous_trial_repeated = false; % initialize this
for e_idx = 1:n_events
  if e.conds(e_idx)==0 % if NULL then don't do anything, no grids..
    % if it's null, consider the previous trial repeated -- i.e., force yourself to make a new grid
    is_previous_trial_repeated = true;
  else % if not null, the continue
    for cycle_ii = 1:p.time.n_cycles_per_tr
      % get a grid 
      % generate a new grid if one of first n or if it's not a repeat
      if is_previous_trial_repeated % don't allow two repeats in a row (i.e., 3 in a row)
        is_repeated_trial = false;
      else
        is_repeated_trial = rand() < p.probability_of_repeat;
      end
      is_first_n_trials = any(isnan(recent_grids(:))); % defining this functionally rather than w/ counters
      if ~is_repeated_trial || is_first_n_trials
        e.grid_is_repeated(e_idx,cycle_ii) = false;
      
        previously_unused_idxs = setdiff(1:n_els_in_grid,previous_grid_idxs);
        new_idxs = randsample(previously_unused_idxs,p.n_bits_flipped,false);
        old_idxs_that_survived = randsample(previous_grid_idxs,p.grid.n_filled_in-p.n_bits_flipped,false);
        filled_grid_idxs = cat(1,old_idxs_that_survived,new_idxs');
        previous_grid_idxs = filled_grid_idxs;

        % fill them in in the grid
        this_grid = zeros(p.grid.grid_height_in_squares,p.grid.grid_width_in_squares);
        for filled_grid_idx = filled_grid_idxs
          this_grid(filled_grid_idxs) = 1;
        end
      
      else % it's a repeat, so just take the relevant previous grid
        e.grid_is_repeated(e_idx,cycle_ii) = true;
        this_grid = recent_grids(:,:,end);
      end
     
      % assign!
      e.grids{e_idx,cycle_ii} = this_grid;
      this_grid_array = generateGridArray(this_grid,p);
      e.grid_textures{e_idx,cycle_ii} = Screen('MakeTexture', w, this_grid_array);

      % update the recent grids and then add this grid
      for recent_ii = 1:(p.n_back-1)
        recent_grids(:,:,recent_ii+1) = recent_grids(:,:,recent_ii);
      end
      recent_grids(:,:,1) = this_grid;
      
      is_previous_trial_repeated = is_repeated_trial; % update this
    end % cycle_ii
  end
end

% and at last, just take a single grid to generate 
%   the square where the textures should be plotted
[s1, s2, ~] = size(this_grid_array);
dstRect = [0 0 s1 s2] .* p.grid.scaling_of_texture_square;
windowRect = Screen('Rect', w);
[xCenter, yCenter] = RectCenter(windowRect);
p.dstRect = CenterRectOnPointd(dstRect, xCenter, yCenter);

end


function gridarray = generateGridArray(grid_pattern,p)
% retrieve relevant vars
square_edge = p.grid.square_edge; grid_width = p.grid.grid_width;
grid_height_in_squares = p.grid.grid_height_in_squares;
grid_width_in_squares = p.grid.grid_width_in_squares;

assert(all(size(grid_pattern)==[grid_height_in_squares, grid_width_in_squares]));

% have a 3d array: (x,y,rgb)
gridarray = zeros(grid_width_in_squares * (square_edge+grid_width) + grid_width,...
                     grid_height_in_squares * (square_edge+grid_width) + grid_width,...
                     3);
for square_ii = 1:grid_height_in_squares
  for square_jj = 1:grid_width_in_squares    
    % determine the color of the square
    if grid_pattern(square_ii,square_jj) % fill it in if should be
      this_square_color = p.grid.on_square_color;
    else % otherwise set it to BG color
      this_square_color = p.grid.off_square_color;
    end
    
    % assign the color to the grid
		gridarray_ii = 1 + grid_width + (square_edge+grid_width)*(square_ii-1);
    gridarray_jj = 1 + grid_width + (square_edge+grid_width)*(square_jj-1);
    gridarray(gridarray_ii:gridarray_ii+square_edge-1,...
                gridarray_jj:gridarray_jj+square_edge-1,...
                :) = repmat(reshape(this_square_color,[1 1 3]),[square_edge,square_edge,1]);
  end
end

end


function e = preloadAudio(e,w,p)
% psychportaudio likes transposed clips

drawCenteredText(w,'Loading audio...',p.fg_color);
Screen('Flip',w);

% add "clip" field based on "clip_names" field

unique_clip_names = unique(e.clip_names)';
is_empty = cellfun('isempty',unique_clip_names);
unique_clip_names = unique_clip_names(~is_empty); % remove empty strings...

% load all clips in just once
unique_clips = {};
n_samples_of_quiet_for_isi = round(p.time.acq_time*p.freq);
for clip_ii = 1:length(unique_clip_names)
  clip_raw = audioread([p.stim_dir filesep unique_clip_names{clip_ii}]); 
  unique_clips{end+1} = [clip_raw; zeros(n_samples_of_quiet_for_isi,2)]; %#ok<AGROW>
end

% place clips in appropriate place...
n_events = length(e.conds);
for event_idx = 1:n_events
  event_cond = e.conds(event_idx);
  
  if event_cond==0
    if ~strcmp(e.clip_names{event_idx},'')
      error('something''s wrong with assigning the clip names -- fixation isn''t an empty string');
    end
    e.clips{event_idx} = [];
  else
    clip_fname = e.clip_names{event_idx};
    % clip_fpath = [p.stim_dir filesep clip_fname];
    % raw_audio = audioread(clip_fpath);
    is_uniq_clip = ismember(unique_clip_names,clip_fname);
    raw_audio = unique_clips{is_uniq_clip};
    
    if e.is_quieter_clip(event_idx)
      raw_audio = 10^(-7/20) * raw_audio; % QUIETER BY 7 DB
    end
    e.clips{event_idx} = raw_audio'*10^(0/20); % transpose for PsychPortAudio
  end
end

Screen('FillRect',w,p.bg_color);
Screen('Flip',w);

end


function trigger_time = waitForTrigger(w,p,which_task_to_do)
drawCenteredText(w,sprintf('%s TASK\n\n(Waiting for trigger...)',upper(which_task_to_do)),p.fg_color);
Screen('Flip',w);

[~,~,key_code] = KbCheck(-3);
while ~any(key_code(p.trigger))
  [~,~,key_code] = KbCheck(-3);
end  

Screen('FillRect',w,p.bg_color);
Screen('Flip',w);
trigger_time = GetSecs();
end


function drawCenteredText(win,text,color) % ,x_offset,y_offset)
if nargin < 3; color = []; end
% if nargin < 4; x_offset = 0; end
% if nargin < 5; y_offset = -30; end

% bounding_box = Screen('TextBounds',win,text);
% bounding_box = Screen(bounding_box, Screen('Rect',win));
% x = bounding_box(RectLeft);
% y = bounding_box(RectTop);

DrawFormattedText(win,text,'center','center',color);
end


function e = presentEvent(e,e_idx,w,p,which_task_to_do)

rct = Screen('Rect',w);
cond = e.conds(e_idx);
clip = e.clips{e_idx};

event_start_absolute = GetSecs();
event_offset_absolute = e.absolute_onsets(e_idx) + e.durs(e_idx);

Screen('FillOval',w,p.fg_color,CenterRect([0 0 p.size_dot p.size_dot],rct)); 
Screen('Flip',w);

if cond==0  % i.e., just do nothing, a null period -- no sound and no visual 
  [key_name_from_whole_event,rt_from_whole_event] = waitLookingForKeys(p.response_keys,event_offset_absolute,event_start_absolute,p);
else
  % initialize these for the auditory events
  key_name_from_whole_event = nan;
  rt_from_whole_event = nan;
  whole_event_onset = cputime();

  % 0. Start playing sound 
  % note that if cond==0, there's no sound to play
  if cond~=0
    MyPsychPortAudio('FillBuffer',p.pa_handle,clip);
    PsychPortAudio('Start',p.pa_handle,1,0,1);
  end
 
  for cycle_ii = 1:p.time.n_cycles_per_tr
    this_texture = e.grid_textures{e_idx,cycle_ii};
    visual_trial_key = nan;
    visual_trial_rt = nan;

    % 0. assign the relevant timings
    this_trial_absolute_onset = e.absolute_onsets(e_idx) + (cycle_ii-1) * p.time.vis_trial_length;
    vis_event_offset_absolute = this_trial_absolute_onset + p.time.vis_stim_length;
    this_trial_absolute_offset = vis_event_offset_absolute + p.time.vis_isi;

    % 1. and put up the grid 
    % last two: rotation and interpolation mode (0:=nearest neighbor)
    Screen('DrawTextures', w, this_texture, [], p.dstRect, 0, 0); 
    Screen('Flip', w);
  
    % 2. Wait, with grid on screen
    [this_key_name,this_rt] = waitLookingForKeys(p.response_keys,...
                                                   vis_event_offset_absolute,...
                                                   this_trial_absolute_onset,...
                                                   p);
    if ~isnan(this_key_name)
      visual_trial_key = this_key_name; visual_trial_rt = this_rt;
      if isnan(key_name_from_whole_event) % record only the first instance; if not first then don't record
        key_name_from_whole_event = this_key_name; rt_from_whole_event = cputime() - whole_event_onset;
      end
    end

    % 2. Record any key press; and figure out if it was right; and color screen accordingly
    %   (if you've gotten here, time has expired)
    has_hit_key = ~isnan(visual_trial_key);
    if has_hit_key
      e.key_pressed_visual(e_idx,cycle_ii) = visual_trial_key;
      e.rt_visual(e_idx,cycle_ii) = visual_trial_rt;
    end  
    if e.grid_is_repeated(e_idx,cycle_ii)==has_hit_key
      is_trial_correct = true;
    else
      is_trial_correct = false;
    end
    
    % 2b. Color screen accordingly
    if strcmp(which_task_to_do,'visual') % only giving feedback for the visual task
      if is_trial_correct
        this_dot_color = p.correct_color;
      else
        this_dot_color = p.incorrect_color;
      end
    else % audio task: randomly color the dot
      if rand() < p.audio_task_prop_visual_randomly_pseudocorrect % pseudo correct
        this_dot_color = p.correct_color;
      else
        this_dot_color = p.incorrect_color;
      end
    end
    Screen('FillOval',w,this_dot_color,CenterRect([0 0 p.size_dot p.size_dot],rct)); 
    Screen('Flip',w);

    % 3. Wait, with grid off screen 
    [this_key_name,this_rt] = waitLookingForKeys(p.response_keys,...
                                                   this_trial_absolute_offset,...
                                                   vis_event_offset_absolute,...
                                                   p);
    if isnan(key_name_from_whole_event) && ~isnan(this_key_name) % if haven't hit a key already and did hit a key now
      key_name_from_whole_event = this_key_name;
      rt_from_whole_event = cputime() - whole_event_onset; 
    end
    
  end % cycle_ii

  % Turn off audio
  if cond~=0 % 0 ==> Null condition and so no sound
    PsychPortAudio('Stop',p.pa_handle);
  end
end

% store the result for the auditory task
if strcmp(which_task_to_do,'auditory')
  e.key_pressed_auditory(e_idx) = key_name_from_whole_event;
  e.rt_auditory(e_idx) = rt_from_whole_event;
end

end


function [key_name, rt] = waitLookingForKeys(keys,time_out,start,p,varargin)
% initialize
rt = nan;
key_name = nan;
if ismember('break-on-keystroke',varargin)
  do_break_on_stroke = true;
else
  do_break_on_stroke = false;
end

while GetSecs() < time_out
  [~, secs, key_code] = KbCheck(-3);
  if any(key_code(keys)) 
    rt = secs - start;
    key = find(key_code~=0,1);
    key_name = find(keys==key);
    if do_break_on_stroke
      break
    end
  end
  if any(key_code(p.quit_keys))
    doEndOfRunHouseKeeping(p,e);
    sca();
    error('User quit out with one of quit key codes');
  end
end
end


function doEndOfRunHouseKeeping(p,e) %#ok<INUSD>
fclose(p.dataFid);
PsychPortAudio('Close',p.pa_handle);

Screen('CloseAll'); % close all textures
% ListenChar(0);

mat_fpath = [p.data_fpath(1:end-3) 'mat']; % don't need to check because already checked for data file
if exist('e','var')
  save(mat_fpath,'e','p');
else
  save(mat_fpath,'p');
end

ShowCursor();
sca();
end
