function [grids,targets] = generate_grids(trial_dur,ioi)
% Create a set of grids for a trial.  In the trial, subjects will be
% detecting grid patterns that are identical to a pattern presented before
% the previous one (see Kell & McDermott, 2019).
% Inputs:
% - trial_dur = trial duration (in s)
% - ioi = inter-onset interval between grids (in ms)
% Outputs:
% - grids = set of grids (nrows x ncols x #grids), where each grid contains 1 if
%     it's filled and 0 otherwise
% - targets = specifies which grids are targets, based on a two-back task
% - params = other parameters relating to the grid display
% Nate Zuk (2019), based on Alex Kell's code

% Grids are ncolsxnrows, where each square is either filled or not filled
% In each grid, nsqon squares are on
ncols = 3; % # columns in the grid
nrows = 2; % # rows in the grid
nsqon = 3; % # of squares on at a time

if ~isempty(varargin)
    for n = 2:2:length(varargin)
        eval([varargin{n-1} '=varargin{n};']);
    end
end

%% Determine the number of grids
nsquares = ncols*nrows;
ngrids = ceil(trial_dur/(ioi/1000))+1; % the number of grids to present

%% Randomly generate grid patterns
squares_on = zeros(nsquares,ngrids);
for ii = 1:ngrids
    sq_idx = randperm(nsquares,nsqon); % randomly select nsqon squares to turn on
    squares_on(sq_idx,ii) = 1; % identify the indexes of squares that should be on
end

%% Determine targets
% Use the grid as a binary representation, and convert to decimal
dec_val = NaN(ngrids,1);
for ii = 1:ngrids
    dec_val(ii) = sum(squares_on(:,ii).*2.^(0:nsquares-1)');
end
% Use the decimal representation to check if there are repeats
targets = [false; false; dec_val(3:end)==dec_val(1:end-2)];

%% Create the grids
grids = false(nrows,ncols,ngrids);
for ii = 1:ngrids
    % convert to logical
    bin_val = logical(squares_on(:,ii));
    % wrap to the grid array
    grids(:,:,ii) = reshape(bin_val,[nrows, ncols]);
end