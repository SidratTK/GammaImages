% independent code to get best Gabor fit for an image at the RF of a chosen
% electrode.
% input: imgRGB = scalar image number from stimulus list 
%                 or array Res x Res x 3 (RGB image).
%                 or array Res x Res. Value (0-1) of hsv h=0, s=1.
%                 or [] = default 1st stimulus image from list
%        degppix = visual degree per pixel of image
%                  if [], loads preprocessed val from stimlus file
%        subject = monkey whose electrode location to use, default alpaH
%        elec = electrode number whose RF location to be used, default 85 ecog
%               or [a e r], 3 element vector with Azi, Ele, Radius in Pixels
% STK 031220

function [imgis,filt1,mP,imgResp] =filterSingleImage(subject,elec,imgIn,degppix,showflag)

% these paths will be used
rootPath = gammaModelPath_st();
dataDir  = fullfile(rootPath, 'Data');

% the next parameters have been decided after checking all of them in rankImgsGammaGui.m
% these are parameters used in the application of filters
widthC        = 3;  % pixel width of circle drawn to mark RF 
widthPatch    = 3;  % display x deg around location 
methodFlag    = '';  % '', construct filters of specified size
getBasicFlag  = 1;   % 1= get filters of same bandwidth from toolbox 
gratingFlag   = 0;   % 0= gabor, 1= grating filter. to be used if method = ''
normFlag      = 1;   % 0= filters in range [-0.5 0.5], 1= filter response normalised
ntypeFlag     = 1;   % 0= norm by best grating response, 1= norm by filtnorm
chooseSize    = [];  % size to check products in SFOri. applied in case getBasicFlag=0. []= use all sizes. 
combinePh     = 1;   % check product after combining phases or use separate phases
nflagOut      = 1;   % scale output response value by best grating response
% getMaxProd    = 1;   % use max filtered product instead of gaussian fit to params
% getPhFlag     = 1;   % also compute the phase
% getSzFlag     = 1;   % also compute best size
% normSzFlag    = 1;   % use normalised filters to compute best size
% isscreenshot  = 1;  % use the screenshot of stimuli

% the following inputs need to be given
if ~exist('subject','var') || isempty(subject), subject = 'alpaH'; end
if ~exist('elec','var') || isempty(elec),       elec = 85;         end
if ~exist('imgIn','var')|| isempty(imgIn),      imgIn= 1;          end
if ~exist('degppix','var'),                     degppix=[];        end
if ~exist('showflag','var'),                    showflag=1;        end
% 1. about the stimulus
if size(imgIn,3)==3      % 3D RGB image provided
    imgRGB = imgIn;
    imghsv = rgb2hsv(imgRGB); 
    imgVal = imghsv(:,:,3)-0.5;  % bring to [-0.5 0.5]
end
if ~isscalar(imgIn) && size(imgIn,3)==1 % 2D Value image provided
    imghsv = cat(3,zeros(size(imgIn)), ones(size(imgIn)), imgIn); % HSV at h=0,s=1, BW
    imgRGB = hsv2rgb(imghsv);
    imgVal = imgIn-0.5;          % bring to [-0.5 0.5]
end
if isscalar(imgIn),              % load processed images?
    chooseim  = imgIn;
    fileName1 = fullfile(dataDir,'allStimuli_noDS.mat');
    load(fileName1, 'allVals','allHues', 'allSats','ImageLabels', 'imageIndsCat', 'Categories');
    imghsv = cat(3,allHues(:,:,chooseim), allSats(:,:,chooseim), allVals(:,:,chooseim));
    imgRGB = hsv2rgb(imghsv);
    imgVal = allVals(:,:,chooseim)-0.5;  % bring to [-0.5 0.5]
    clear allVals allHues allSats;       % for memory saving
end
if isempty(degppix)
    fileName1 = fullfile(dataDir,'allStimuli_noDS.mat');
    load(fileName1, 'degppix');
end
res = size(imgRGB,1);

% 2. Load RF details & get location in pix
[rfStatsDeg,~,LFPElectrodeList,EcogElectrodeList] = getRFdetails(subject,dataDir);
rfStatsPix = changeRFDegtoPix(rfStatsDeg,degppix,res);    % convert to pixels
rfStats    = cell2mat(rfStatsPix);    % unwrap from cell

if length(elec)~=3    % if location in pix not provided
   elec = elec(1);    % extra check. pick electrode number
   a = rfStats(elec).meanAzi;
   e = rfStats(elec).meanEle;
   r = sqrt(0.5*(rfStats(elec).rfSizeAzi^2 + rfStats(elec).rfSizeEle^2));  % RF rad in pixels. see DubeyRay 2019JNS p4301
else
   a = elec(1); e= elec(2); r = elec(3);
end
paramsRF(1) = a;         % RF center in pixels
paramsRF(2) = e;
paramsRF(3) = r;         % make a circle at RF size


% 3. get the best Gabor for the image patch
[mP, imgResp, spatFilter] = getBestGaborFit(imgVal,a,e,res,degppix,chooseSize,methodFlag,normFlag,ntypeFlag,gratingFlag,getBasicFlag,combinePh,nflagOut);
filtis = spatFilter{1}{1}+0.5;    % unwrap & get bw 0-1

% 4. Mark RF on the image & selected filter 
imgis = insertShape(imgVal+0.5,'Circle', paramsRF(1:3), 'LineWidth', widthC, 'Color','r');                
filt1 = insertShape(filtis,'Circle', paramsRF(1:3), 'LineWidth', widthC, 'Color','r');                

% 5. display cropped section around electrode RF 
if showflag
    figure; 
    implots(1) = subplot(211); implots(2) = subplot(212); 
    makeImPltLims(a,e,widthPatch,degppix,implots(1),imgis);
    makeImPltLims(a,e,widthPatch,degppix,implots(2),filt1);
end

end 

function [mP,imgResp,spatFilterOut] = getBestGaborFit(imgVal,rfcol,rfrow,res,degppix,chooseSize,methodFlag,normFlag,ntypeFlag,gratingFlag,getBasicFlag,combinePh,nflagOut)

DegIm = res*degppix;                   % expanse of image

% for making filters
SFs    = [0.5 1 2 4 8];                % cycles per vis. degree
Sizes  = [0.1 0.2 0.4 0.8 1.6 3.2];    % visual degrees. sigma or radius/3. vals used in Size ori
Thetas = deg2rad([0:22.5:157.5]);
Phases = [0, pi/2];

% make required type of filters
[spatFilter, filterWidth, ~,~,~,optimGratResp,normoffilt] = makeGratingFilter(SFs, Thetas, Sizes, Phases, gratingFlag, res, degppix, 0,getBasicFlag);

% move them to location
[spatFilterFull,~] = moveGratingFilter(rfcol,rfrow,spatFilter,filterWidth);
         
% get the products
imgPrd{1} = filterRFimage(spatFilterFull,imgVal,rfcol,rfrow);
filteredProd = imgPrd;       % save in a cell so its compatible with next fnxns

% normalise products if required
if normFlag && ~strcmp(methodFlag,'0') 
    if ntypeFlag==1
        filteredProd = cellfun(@(x) x./normoffilt,filteredProd,'un',0);
    else
        filteredProd = cellfun(@(x) x./optimGratResp,filteredProd,'un',0);
    end
end

% quad phase sum of 2 filter responses.
if combinePh && ~strcmp(methodFlag,'0') 
    filteredProd = cellfun(@(x) sqrt(x(:,:,:,1).^2 + x(:,:,:,2).^2), filteredProd,'un',0);
end

% use the SFOri plane at chosen size if selected
if getBasicFlag
    useSizes=0;       % because filter sizes depend on SF, 0 is used in checkMaxProd for getBasic
else
    useSizes = Sizes;
    if ~isempty(chooseSize)
        useSizeInd  = getClosestInd(Sizes,chooseSize);
        useSizes    = Sizes(useSizeInd);
        filteredProd= cellfun(@(x) x(useSizeInd,:,:,:), filteredProd,'un',0);
    end  
end
    
% find max product and get Sf, ori
[mP, ~] = checkMaxProd(filteredProd,filterWidth,SFs,Thetas,useSizes,Phases,degppix);
sP = nan(size(mP));   % mP= meanParameters (or best fit), sP= sigma but not used in present case. 

% calculate phase
mP = getDFTPh(imgVal,mP,degppix,DegIm,rfcol,rfrow,getBasicFlag,imgPrd,useSizes,SFs,Thetas);        
   
% calculate size
[mP,~] = getSizeChecked2(imgVal,degppix,rfcol,rfrow,Sizes,mP,sP,1,ntypeFlag);

% make output filter and get normalised response
if isnan(mP(4,1)),    ph = 0;    else    ph = mP(4,1);    end
[spatFilter,filterRad,~,~,~,optResp,~] = makeGratingFilter(mP(2), mP(3), mP(1), ph, gratingFlag, res, degppix, 0, 0);
[spatFilterOut,~] = moveGratingFilter(rfcol,rfrow, spatFilter, filterRad);

% apply to image:
imgResp = filterRFimage(spatFilterOut,imgVal,rfcol,rfrow);
if nflagOut
    imgResp = imgResp./optResp; % always wrt optimum grat resp, so bw 0-1, signifies contrast
end 

end

function [MaxP, AmpP] = checkMaxProd(filteredProd,filterWidth,SFs,Thetas,Sizes,Phases,degppix)

numStimuli = length(filteredProd);

% for each image, find the max product with the filter 
for im = 1:numStimuli
    temp = abs(filteredProd{im});     % magnitude of product
    if size(temp,4)>1                 % if both phases product available
        [X,Y,Z,W] = ndgrid(Sizes,SFs,Thetas,Phases);
    else                              % if products combined
        [X,Y,Z]   = ndgrid(Sizes,SFs,Thetas);
    end
    yD   = temp(:);
    [AmpP(im), ind] = max(yD); % amp is the amplitude.
    MaxP(1,im) = X(ind);       % size  
    MaxP(2,im) = Y(ind);       % sf
    MaxP(3,im) = Z(ind);       % theta
    if size(temp,4)==2 
        MaxP(4,im) = W(ind);
    else
        MaxP(4,im) = nan;
    end
    if MaxP(1,im)==0 && ~isempty(filterWidth)  % will happen when getBasic=1 % read from filterWidth
        indis = getClosestInd(SFs,MaxP(2,im));
        MaxP(1,im) = (filterWidth(1,indis)*degppix)/3;  % get sigma in degrees
    end
        
end 
    
end

function [mP]= getDFTPh(allVals,mP,degppix,degIm,azi,ele,getBasic,imgProd,Sizes,SFs,Thetas)
% allVals has images
% mP has all other features

defaultsize = 1;  % if none provided, use this sigma to Gaussian mask image section we want. 
res = size(allVals,1);

% ready for a symmetrical Gaussian envelope at location
[columnsInImage, rowsInImage] = meshgrid(1:res, 1:res);

for im = 1:size(allVals,3)
   if isnan(mP(1,im)) || isempty(mP(1,im)) || mP(1,im)==0  % should not arise ideally
      mP(1,im) = defaultsize;
   end  
   sd = mP(1,im)/degppix;
   gaus = exp( ((columnsInImage-azi).^2+(rowsInImage-ele).^2)/-(2*(sd)^2) );   % ranges from 0-1
    
   % ph needs to be corrected. first, the phase detected for a 0 phase grating centered at RF?
   [spatGab,isRad,~,~,~,~] = makeGratingFilter(mP(2,im), mP(3,im), mP(1,im), 0, 0, res, degppix, 0,getBasic); % large size   
   [spatGabOut]= moveGratingFilter(azi, ele, spatGab, isRad);
   ph0grating  = spatGabOut{1}{1};
   [ph0] = getPhaseSelected(ph0grating, degppix, degIm, mP(2,im), mP(3,im));
   
   useVal = allVals(:,:,im).*gaus;  % smooth fade out the images using gaus
   [ph]   = getPhaseSelected(useVal, degppix, degIm, mP(2,im), mP(3,im)); % get DFT phase at location 
  
   mP(4,im) = ph-ph0;
   
end

% % try doing same using the 0 and pi/2 phase products available
% % but atan has range -pi/2 to pi/2. It mistakes the 2nd and 3rd quadrant
% for im = 1:size(allVals,3)
%     useProd =  imgProd{im};
%     indS = getClosestInd(Sizes,mP(1,im));
%     indF = getClosestInd(SFs,mP(2,im));
%     indT = getClosestInd(Thetas,mP(3,im));
%     
%     useProd = squeeze(useProd(indS,indF,indT,:));
%     phout1(im) = atan(useProd(2)/useProd(1));
% end
% mP(4,:) = phout1;

end

function [mPOut,sPOut,prodOut] = getSizeChecked2(allVals,degppix,azi,ele,Sizes,mP,sP,normSz,ntypeFlag)
mPOut = mP;
sPOut = sP;

res = size(allVals,1);
for im = 1:size(allVals,3)
    SF = mP(2,im);
    Or = mP(3,im);
    Ph = mP(4,im); if isnan(Ph), Ph=0; end
    % make filters
    [FiltsAre,radAre,~,~,~,optGratResp,filtNorm] = makeGratingFilter(SF,Or,Sizes,Ph,0,res,degppix,0,0);
    FiltsAre = moveGratingFilter(azi, ele, FiltsAre, radAre);
    tempPrd  = filterRFimage(FiltsAre,allVals(:,:,im),azi, ele);
    if normSz
        if ntypeFlag, tempPrd = tempPrd./filtNorm;
        else          tempPrd = tempPrd./optGratResp;
        end 
    end
    
    tempPrd = round(tempPrd,2);  %  approximation at 1/100th level
    indSi      = find(tempPrd == max(tempPrd));
    mPOut(1,im)= max(Sizes(indSi));   % choose largest size if more available 
    sPOut(1,im)= nan;
    prodOut(1,im)= max(tempPrd);
end 

end
