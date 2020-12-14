% code to move a given grating filter to desired location.
% 021020; used floor instead of round to round the desired rows and columns
% since filterPixels-1 steps are taken towards the right, so floor keeps
% the center closest to the RF (which can be at 330.4 pixels eg) across diff images
% 141020: if i take floor, it should be filterPixels-1 to the left, not
% right. so that the same no of pixels are chosen on eithe side of the RF. 
% 141020: updated to send in resout explicitly & o/p size not same as temp


function [spatFilterOut,spatFilterCut]= moveGratingFilter(rfcol, rfrow, spatFilter, filterRadius, resout)
% rfcol and rfrow are wrt to output filter.

numSizes = size(spatFilter,1);
numSFs = size(spatFilter,2);
numThetas = size(spatFilter,3);
numPhases = length(spatFilter{1});

spatFilterOut = cell(size(spatFilter));
spatFilterCut = cell(size(spatFilter));
% displace the grating from center to the RF center.
for s = 1:numSizes
    for f = 1:numSFs
        filterPixels = filterRadius(s,f);
        for t = 1:numThetas
            for ph = 1:numPhases
                temp = spatFilter{s,f,t}{ph};
                res  = size(temp,1);       % resolution
                centerpix = (1+res)/2;
                if ~exist('resout','var'),
                    resout = size(temp,1);
                end
                spatFilterOut{s,f,t}{ph} =  zeros(resout,resout);
                
                pixSourceRows = floor([centerpix-(filterPixels-1) : centerpix+filterPixels]);
                pixSourceCols = floor([centerpix-(filterPixels-1) : centerpix+filterPixels]);
                
                spatFilterCut{s,f,t}{ph} = temp(pixSourceRows,pixSourceCols);
                
                pixDestRows   = floor([rfrow-(filterPixels-1) : rfrow+filterPixels]);
                pixDestCols   = floor([rfcol-(filterPixels-1) : rfcol+filterPixels]);
                % if filter going out of image edges, truncate it in that case
                pixSourceRows = pixSourceRows(pixDestRows>0 & pixDestRows<=resout);
                pixSourceCols = pixSourceCols(pixDestCols>0 & pixDestCols<=resout);
                pixDestRows   = pixDestRows(pixDestRows>0 & pixDestRows<=resout);
                pixDestCols   = pixDestCols(pixDestCols>0 & pixDestCols<=resout);
                
                spatFilterOut{s,f,t}{ph}(pixDestRows,pixDestCols) = temp(pixSourceRows,pixSourceCols);
            end
        end 
    end
end


end
