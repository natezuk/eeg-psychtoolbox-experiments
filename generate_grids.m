function [grids,targets] = generate_grids(trial_dur,ioi)
% Create a set of grids for a trial.  In the trial, subjects will be
% detecting grid patterns that are identical to a pattern presented before
% the previous one (see Kell & McDermott, 2019).
% Inputs:
% - trial_dur = trial duration (in s)
% - ioi = inter-onset interval between grids (in ms)
% Outputs:
% - grids = set of 3x3 grids (3x3x#grids), where each grid contains 1 if
%     it's filled and 0 otherwise
% - targets = specifies which grids are targets, based on a two-back task
% - params = other parameters relating to the grid display
% Nate Zuk (2019), based on Alex Kell's code

% Grids are 3x3, where each square is either filled or not filled
ncols = 3; % # columns in the grid
nrows = 2; % # rows in the grid
nsquares = ncols*nrows;

%% Determine the number of grids
ngrids = ceil(trial_dur/(ioi/1000))+1; % the number of grids to present

%% Randomly generate grid patterns
% Start by generating random numbers between 0 and 2^8-1, since each grid can
% be represented as one of these numbers (binary on-off for each grid
% square)
grid_vals = randi(2^(nsquares-1)-1,ngrids,1);

%% Determine targets
targets = [false; false; grid_vals(3:end)==grid_vals(1:end-2)];

%% Create the grids
grids = false(nrows,ncols,ngrids);
for ii = 1:ngrids
    % convert to a binary number, array of 1 and 0
    bin_val = dec_to_bin_array(grid_vals(ii),nsquares);
    % convert to logical
    bin_val = logical(bin_val);
    % wrap to the grid array
    grids(:,:,ii) = reshape(bin_val,[nrows, ncols]);
end