% making the filter bank
% updated 290420 to 1. remove dependency on knkutils.
%                   2. radius = 3*sigma
%                   3. gratings made correspond to those shown in Knot
% updated 080520 to make Gabor filters
% updated 180520 to add normalisation of filter by response to optimum full size grating 
% updated 070620 to make filter using KNK toolbox as well using useTBflag
% updated 180820 to remove optimRespOut output. It was causing problem in toolbox filter
% updated 180820 to round radXY in make2dgabor
% updated 310820 deleted extra fnxn getResfromImg. it was never used.
% updated 030920 deleted threshold criteria in make2dGabor, which made gaussian values below thresh 0. 
% As the 3sigma=radius criteria is used to limit the Gabor anyway, These steps are redundant.
% 060920: observed that the filters from here (for ph=0) are 90 degrees shifted
% compared to the knot stimuli. but this is no problem since it is always 90. 
% make2dGratingFull removed & saved in functions required.
% 220920: filt_prop.spacings increased so that memory constraints down the fnxns are not a problem
% 240920: in makeGaborfromTB, do normalisation as in makeSpatialGrating fnxn, not via filt_prop.scaling parameter.
% 131020: outputs the filtnormOut whihc is the self norm of the filter.
% Does not normalise filters before hand. 
% 141020: call moveGratingFilter in makeGaborfromTB 

% STK 

function [spatFilter,filterRadius,centerpix,res,degppix,optimRespOut,filtnormOut] = makeGratingFilter(SFs, Thetas, Sizes, Phases, gratingFlag, res, degppix, showFlag,useTBflag)
if ~exist('gratingFlag','var'),   gratingFlag = 0;   end % gabor or grating
if ~exist('isscreenshot','var'),  isscreenshot= 1;   end % edited screenshots or full images
if ~exist('showFlag','var'),      showFlag    = 1;   end
if ~exist('useTBflag','var'),     useTBflag   = 0;   end % use the tool box codes of Knk?
if ~exist('res','var'),           res       = 540;   end % pixel dimension of sq image
if ~exist('degppix','var'),       degppix= 0.0463;   end % visual degree per pixel

if ~exist('SFs','var'),     SFs    = [0.5 1 2 4 8];             end    % cycles per vis. degree
if ~exist('Thetas','var'),  Thetas = deg2rad([0:22.5:157.5]);   end    % angular degrees. vals used in Size ori
if ~exist('Sizes','var'),   Sizes  = [0.1 0.2 0.4 0.8 1.6 3.2]; end    % visual degrees. sigma or radius/3. vals used in Size ori
if ~exist('Phases','var'),  Phases = [0, pi/2];                 end


numSFs = length(SFs);
numThetas= length(Thetas);
numSizes = length(Sizes);
numPhases= length(Phases);
DegIm    = degppix*res;

if ~useTBflag
    [spatFilter, filterRadius, centerpix, optimRespOut,filtnormOut] = makeSpatialGrating(SFs,Thetas,Sizes,Phases,degppix,res,gratingFlag); 
    % spatFilter 3D cell. size x sf x theta. Each cell has 2 quadrature phase Gabors/ gratings
    filterRadius = repmat(filterRadius(:),1,length(SFs));  % in pixels
    
else   % use tool box of KNK paper. Size parameter irrelevant
    [spatFilter, filterRadius, centerpix, optimRespOut,filtnormOut] = makeGaborfromTB(SFs,Thetas,Phases,DegIm,res,degppix);
    filterRadius = repmat(filterRadius(:)',length(Sizes),1); 
end

if showFlag   % display Gratings/Gabors. 
    figure; colormap gray;
    plotsA = getPlotHandles(numSFs, size(spatFilter,1), [0.1 0.1 0.8 0.8], 0.04,0.05);
    for f = 1:numSFs 
        for s = 1:size(spatFilter,1)
            t = randperm(numThetas,1); % pick an ori
            ph= randperm(2,1); 
            gabis = spatFilter{s,f,t}{ph};
            % plot here
            axes(plotsA(f,s));
            imagesc(gabis);
            caxis([-0.5 0.5]);
            xlabel(plotsA(f,s), [' Ang ',num2str(rad2deg(Thetas(t)))],'fontWeight','bold');    
            if f==1
                title(plotsA(f,s),['Size ',num2str(Sizes(s))],'fontWeight','bold');    
            end 
        
        end
        ylabel(plotsA(f,1), ['SF ',num2str(SFs(f)),'cpd'],'fontWeight','bold');
    end
end 

end

function [spatFilter,filterRad,centerpix,optimRespOut,filtnormOut] = makeSpatialGrating(SFs,Thetas,Sizes,Phases,degppix,res,gratingFlag)

% SFs : the spatial frequencies of grating
% Thetas: angle in radians, 0 = horizontal grating 
% Sizes : is the sigma of the Gaussian or the radius of the aperture.
% degppix: visual degrees per pixel 
% gratingFlag: 0 = make Gabor by Gaussian envelope. 1 = use discrete
% aperture

if ~exist('gratingFlag','var'), gratingFlag= 0; end
if gratingFlag
    Sizes = Sizes*3;    % make radius of grating from sigma of gabor. r/s=3.
end

numSFs = length(SFs);
numThetas= length(Thetas);
numSizes = length(Sizes);
numPhases= length(Phases);

spatFilter= cell(numSizes,numSFs,numThetas);

for s = 1:numSizes
    for f = 1:numSFs
        for t = 1:numThetas
            for ph = 1:numPhases
                size = Sizes(s); sf = SFs(f); ang = Thetas(t); phase = Phases(ph);
                
                if gratingFlag  % make Grating  
                    [gratingIs,~,radXY,centerpix] = make2dGrating(res,degppix,size,sf,ang,phase); 
                else            % make Gabor
                    [gratingIs,~,radXY,centerpix] = make2dGabor(res,degppix,size,sf,ang,phase);   
                end
                % it is made in range -1 to 1. Bring it to -0.5 to 0.5
                gratingIs = gratingIs/2;
                 
                % there are 2 ways I can normalise the filter if needed: 
                % a) by the response to the optimum grating as done by Kedrick Kay toolbox 
                optimGrat = (make2dGratingFull(res,degppix,sf,ang,phase))/2;
                optimRespOut(s,f,t,ph) = sum(sum(optimGrat.* gratingIs)); % similar to knk paper applymultiscalegabors
                
                % b) by the norm of the filter itself
                filtnormOut(s,f,t,ph) = sqrt(sum(sum(gratingIs.^2)));
                
                % not normalising here, will do later when reading products
%                 if normFlag       % normalise the filter by response to optimal grating
%                     gratingIs = gratingIs/optimRespOut(s,f,t,ph);
%                 end

                spatFilter{s,f,t}{ph} = gratingIs;
            end      
        end 
    end 
    filterRad(s) = radXY; 
end

end

function [gaborIs,gratingFull,radXY,r] = make2dGabor(res,degppix,size,sf,ang,phase)

sd = size/degppix;   % convert sigma of the gabor from degree to pixel
thresh = 0.01;       % Gaussian values below this made 0. optional. commented usage

% identify center:
r = (1+res)/2;     % central row
c = (1+res)/2;     % central column

% get a grating 
[gratingFull] = make2dGratingFull(res,degppix,sf,ang,phase);

% make a symmetrical Gaussian envelope at center
[columnsInImage, rowsInImage] = meshgrid(1:res, 1:res);
gaus = exp( ((columnsInImage-c).^2+(rowsInImage-r).^2)/-(2*sd^2) );   % ranges from 0-1
% gaus(gaus<thresh) = 0;    % commented threshold criteria

% apply mask:
gaborIs  = gratingFull.*(gaus);

% further crop the Gabor at 3*sigma by using circular mask
radXY = round(3*sd);  % radius taken as 3 times sigma. pixels 
ellipsePixels = ((rowsInImage-r).^2 ./radXY^2 + (columnsInImage-c).^2 ./radXY^2) <=1 ;  % 2d logical array 
ellipsePixels = cast(ellipsePixels,'like',gaborIs);
% apply mask:
gaborIs  = gaborIs.*(ellipsePixels);

end

function [gratingIs,gratingFull,radXY,r] = make2dGrating(res,degppix,size,sf,ang,phase)

% get a grating first
[gratingFull] = make2dGratingFull(res,degppix,sf,ang,phase);

% identify center:
r = (1+res)/2;     % central row
c = (1+res)/2;     % central column

% Now make a circular mask:
radXY = round(size/degppix);  % radius in pixels 
[columnsInImage, rowsInImage] = meshgrid(1:res, 1:res);
ellipsePixels = ((rowsInImage-r).^2 ./radXY^2 + (columnsInImage-c).^2 ./radXY^2) <=1 ;  % 2d logical array 
ellipsePixels = cast(ellipsePixels,'like',gratingFull);

% apply mask:
gratingIs  = gratingFull.*(ellipsePixels);

end

function [spatFilter,filterRadius,centerpix,optimRespOut,filtnormOut] = makeGaborfromTB(SFs,Thetas,Phases,DegIm,res,degppix)

% put thetas & phases in cell for modified applymultiscalegaborfilters function
PhasesIn = {Phases};  
ThetasIn = {Thetas};  
numThetas = length(Thetas);
numPhases = length(Phases);
filt_prop.spatFreq    = SFs;
filt_prop.cycles      = SFs* DegIm;   % number of cycles per image field of view
filt_prop.bandwidth   = -1;           % the spatial frequency bandwidth of the filters is 1 octave
filt_prop.orientations= ThetasIn;     % filters occur at 8 orientations
filt_prop.phases      = PhasesIn;     % filters occur at 2 phases (between 0 and pi)
filt_prop.thres       = 0.01;         % the Gaussian envelopes are thresholded at .01
filt_prop.mode        = 0;            % the dot-product between each filter and each stimulus is computed
filt_prop.spacings    = 25;           % large spacing bw Gabors because i only need the filter here, not products
filt_prop.scaling     = [0 0.5];      % filts bw [-0.5 0.5]. % 2=>filters scaled to have an equivalent Michelson contrast of 1

dummyImg = zeros(res,res);         % dummy for next function
    
    % using applymultiscalegaborfilters so that gabors get normalised and cropped as well
[~,gbrs,~,sds,~,~,~] = applymultiscalegaborfilters(reshape(dummyImg(:,:,1),res*res,[])', ...
    filt_prop.cycles,filt_prop.bandwidth,filt_prop.spacings,filt_prop.orientations,...
    filt_prop.phases,filt_prop.thres,filt_prop.scaling,filt_prop.mode);

centerpix    = (1+res)/2;
filterRadius = cell2mat(cellfun(@(x) 0.5*size(x,1), gbrs(:,1,1),'un',0));  % same for all thetas and phases 
for s = 1  % size depends on SF
    for f = 1:length(SFs)
        for t = 1:numThetas
            for ph = 1:numPhases
                tempFilter{1,1,1}{1} = gbrs{f,t,ph};
                [tempFilterfull,~] = moveGratingFilter(centerpix,centerpix,tempFilter,filterRadius(f),res);
                gratingIs = tempFilterfull{1}{1};
                
                % a) can normalise by response to optimum grating
                optimGrat = (make2dGratingFull(res,degppix,SFs(f),Thetas(t),Phases(ph)))/2;
                optimRespOut(s,f,t,ph) = sum(sum(optimGrat.* gratingIs)); % similar to knk paper applymultiscalegabors
                
                % b) or by the norm of the filter itself
                filtnormOut(s,f,t,ph) = sqrt(sum(sum(gratingIs.^2)));
                
                % not normalising here, will do later when reading products
%                 if normFlag       % normalise the filter by response to optimal grating
%                     gratingIs = gratingIs/optimResp;
%                 end
                spatFilter{s,f,t}{ph} = gratingIs;
                
            end 
        end
    end
end  

% I have checked the above filters (that are moved) with the ones that are made by
% applymultiscalegaborfilters. They are same if indices are same. 
 
end

