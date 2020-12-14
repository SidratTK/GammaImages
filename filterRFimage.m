

function [filteredProd] = filterRFimage(spatFilter,stimis,rfcol,rfrow)
% spatFilterSegment has cells with a 2 D nxn, matrix with spatial filter values
% stimis is 2D mxm, matrix with image values
% n<m
% STK 050420

rfrow= round(rfrow); rfcol=round(rfcol);   % if between two pixels
numSizes = size(spatFilter,1);
numSFs = size(spatFilter,2);
numThetas = size(spatFilter,3);
numPhases = length(spatFilter{1});

for s = 1:numSizes
    for f = 1:numSFs
        for o = 1:numThetas
            for ph = 1:numPhases
                filteris = spatFilter{s,f,o}{ph};
                
                if size(filteris) == size(stimis)     % size matched filter
                    temp = sum(sum(filteris.*stimis));
                    filteredProd(s,f,o,ph) = temp;
                elseif size(filteris) < size(stimis)  % segmented filter
                    filtImage = filter2(filteris,stimis,'same');         % apply convolution based filter.
                    filteredProd(s,f,o,ph) =  filtImage(rfrow,rfcol);    % pick product at RF center 
                elseif size(filteris) > size(stimis) 
                    disp('filter is larger than image. Unpredicted case should not arise. see also moveGratingFilter');
                end   
                
            end 
        end
    end
end

end

