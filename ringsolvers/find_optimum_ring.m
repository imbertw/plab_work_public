function [opt_rOut opt_width opt_FOM] = find_optimum_ring(rOut, width, figure_of_merit)
[row,column] = find(figure_of_merit==max(max(figure_of_merit)));
opt_rOut = rOut(row);
opt_width = width(column);
opt_FOM = figure_of_merit(row,column);