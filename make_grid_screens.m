function grid_screens = make_grid_screens(scnPtr,grids)
% Create a set of screens for Psychtoolbox containing grids, which will be
% presented for a 2-back visual task
% Nate Zuk (2019)

%%% Position a crosshair in the center, one row of squares on the top, and
%%% one row of squares on the bottom

upper_grid_position = 0;
lower_grid_position = 0;
cross_hair_color = 0;
cross_hair_position = 0;
cross_hair_size = 0;
on_square_color = [];
off_square_color = [];
