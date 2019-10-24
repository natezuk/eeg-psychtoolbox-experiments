% Code from Alex Kell's file responsible for generating and displaying
% grids for the visual task

%% Grid setup
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


%% Display the grid and keep track of responses

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

%% Grid design
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