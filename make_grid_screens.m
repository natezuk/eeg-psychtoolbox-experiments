function grid_textures = make_grid_screens(wS,grid_patterns,varargin)
% Create a set of screens for Psychtoolbox containing grids, which will be
% presented for a 2-back visual task
% Nate Zuk (2019)

%%% Position a crosshair in the center, one row of squares on the top, and
%%% one row of squares on the bottom

% upper_grid_position = 0;
% lower_grid_position = 0;
square_width = 9;
square_edge_width = 1;
crosshair_active = 1; % flag to draw the crosshair
crosshair_color = [255 255 255];
lines_color = [255 255 255]; % lines of the grid, in [r g b]
background_color = wS.bckColor;
off_square_color = background_color;
on_square_color = [204 0 204]; % magenta

% Parse varargin (to change default parameters)
if ~isempty(varargin)
    for n = 2:2:length(varargin)
        eval([varargin{n-1} '=varargin{n};']);
    end
end

% Get the height and width of the grid in squares
grid_height = size(grid_patterns,1);
grid_width = size(grid_patterns,2);
grid_textures = cell(size(grid_patterns,3),1);

% Go through each grid
for n = 1:size(grid_patterns,3)
    %% Make grid
    % Setup a 3-D matrix representing the pixels of the texture (x,y,rgb)
    % (courtesy of Alex Kell)
    viewing_corner_pxl = zeros(grid_height,grid_width,2); % defines the top-left-most
        % corner (x,y coordinate) of the square corresponding to the identical index 
        % in the grid
    % Adjust the grid to place the crosshair in the center
    % (crosshair is defined to be slightly smaller than the square width)
    % If there are an even number of squares in the grid, place squares
    % evenly around the crosshair so the crosshair is visible at the center
    if mod(size(grid_pattern,2),2)==0 && crosshair_active % if # columns are even
        viewed_grid_width = (grid_width+1) * (square_edge+square_edge_width) + square_edge_width;
        % identify the corner index within the viewed where a square should be
        % allocated from the grid pattern
        grid_idx = [1:floor((grid_width+1)/2) ceil((grid_width+1)/2)+1:(grid_width+1)]; % skip the center index
    else
        viewed_grid_width = grid_width * (square_width+square_edge_width) + square_edge_width;
        grid_idx = 1:grid_width;
    end
    % compute the top-left x-coordinate for each square
    x_corner_pxls = 1 + square_edge_width + (square_edge+square_edge_width)*grid_idx;
    viewing_corner_pxl(:,:,1) = repmat(x_corner_pxls,grid_height,1); % copy for each square with the same x coord
    if mod(size(grid_pattern,1),2)==0 && crosshair_active % if # rows are even
        viewed_grid_height = (grid_height+1) * (square_edge+square_edge_width) + square_edge_width;
        % get the corner y-coord for each square
        grid_idx = [1:floor((grid_height+1)/2) ceil((grid_height+1)/2)+1:(grid_height+1)]; % skip the center index
    else
        viewed_grid_height = grid_height * (square_edge+square_edge_width) + square_edge_width;
    end
    % compute the top-left y-coordinate for each square
    y_corner_pxls = 1 + square_edge_width + (square_edge+square_edge_width)*grid_idx';
    viewing_corner_pxl(:,:,2) = repmat(y_corner_pxls,1,grid_width);

    % preallocate the viewing grid, and set color to line color
    gridarray = NaN(viewed_grid_width,viewed_grid_height,3);
    for ii=1:3 % for each r,g,b
        gridarray(:,:,ii) = lines_color(ii);
    end

    for square_ii = 1:grid_height
      for square_jj = 1:grid_width   
        % determine the color of the square
        if grid_patterns(square_ii,square_jj) % fill it in if should be
          this_square_color = on_square_color;
        else % otherwise set it to BG color
          this_square_color = off_square_color;
        end

        % assign the color to the grid
        gridarray((0:square_edge-1)+viewing_corner_pxl(:,:,1),...
                  (0:square_edge-1)+viewing_corner_pxl(:,:,2),...
                  :) = repmat(reshape(this_square_color,[1 1 3]),[square_edge,square_edge,1]);
      end
    end


    %% Make the crosshair
    if crosshair_active,
        % Get the center of the grid
        gridcenter = [round(viewed_grid_width/2) round(viewed_grid_height)/2];
        % draw the x span of the crosshair
        cross_xspan = gridcenter(1)+(0:square_edge-2);
        gridarray(cross_xspan,gridcenter(2),:)...
            = repmat(reshape(crosshair_color,[1 1 3]),[square_edge-2,1,1]);
        % draw the y span of the crosshair
        cross_yspan = gridcenter(2)+(0:square_edge-2);
        gridarray(gridcenter(1),cross_yspan,:)...
            = repmat(reshape(crosshair_color,[1 1 3]),[1,square_edge-2,1]);
    end

    %% Draw the grid
    grid_textures{n} = Screen('MakeTexture', wS.ptr, gridarray);
    
end