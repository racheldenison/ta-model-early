function sw = makeSquareWave(space,center,width,height)
%
% squareWave = makeSquareWave(space,center,width,[height])
%
% This function creates a square wave centered at "center", 
% over the values defined by vector "space".  
%
% height, if specified is the height of the square wave.
% Otherwise, it is scaled to unit volume

plot_figure = 0; % if 1, then plot the two kernels.

sw = unifpdf(space,center-width/2,center+width/2); 

if exist('height','var')
  sw(sw > 0) = height;
end

if plot_figure == 1
    figure; plot (space,sw)
end
