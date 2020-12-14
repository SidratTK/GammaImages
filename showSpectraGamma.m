% function to display the freq spectra and
% gamma power in stimulus differnet conditions.
% modified from showSpectraUsed.m 
% STK 031120

function showSpectraGamma()

rootPath = fileparts(which(mfilename));
dataDir  = fullfile(rootPath, 'Data');
folderSourceString = '/Volumes/SIDSANDISK/fss/';
gridType = 'Microelectrode';

useNotchData = 0;   % whether notched before spectra calculation
methodFlag   = 0;   % 0= MT & dB change in gamma, 1= Gaussian gamma fit like Hermes etal.
tBL = [-0.25 0];    % BL period
tST = [0.25 0.5];   % stim period
numTapers = 3; 

subjectNames = {'alpaH','kesariH'};
fBand = [35 70]; % [45 70];   %    % Gamma band 

useoriforsf  = {45 90}; % 2 monkeys, ori for SF tuning
usesfforori  = {4 2};  % sf for ori tuning
useoriforcon = {[] []};  % when [], average over all oris
useoriforsze = {[] []};
% figure related
showForParams = {'Size','SF','Ori','Con'};
stimCols{1}  = parula(6);% size
stimCols{2}  = hsv(6);   % sf
stimCols{3}  = hsv(8);   % ori  
stimCols{4}  = gray(9);  % con
lfpcol = 'm';
ecogcol= 'b';
fTicks   = [0:50:200];
tickLength = [0.015 0]; fontSizeLabel = 13; fontSizeTicks = 9; fontSizeSmall = 10; fontSizeLarge = 11;

% see spectra & gamma tuning 
figH = figure;
figH.Units = 'centimeters';
figH.PaperType  = 'a4';
figH.PaperUnits = 'centimeters';
figH.PaperSize  = [17.5 25]; % figH.PaperSize  = [17.8 25.1];
figH.PaperOrientation = 'Portrait';
figH.PaperPosition = [0 0 figH.PaperSize];
figH.Color = [1 1 1];
figH.Position = [0 0 figH.PaperSize];
figH.PaperUnits = 'normalized';

SizePlots      = getPlotHandles(2,2,[0.07 0.70 0.4 0.24], 0.04, 0.02);
SizePlots(3,:) = getPlotHandles(1,2,[0.07 0.54 0.4 0.10], 0.04, 0.03);
SfPlots        = getPlotHandles(2,2,[0.57 0.70 0.4 0.24], 0.04, 0.02);
SfPlots(3,:)   = getPlotHandles(1,2,[0.57 0.54 0.4 0.10], 0.04, 0.03);
CnPlots        = getPlotHandles(2,2,[0.07 0.21 0.4 0.24], 0.04, 0.02);
CnPlots(3,:)   = getPlotHandles(1,2,[0.07 0.05 0.4 0.10], 0.04, 0.03);
OrPlots        = getPlotHandles(2,2,[0.57 0.21 0.4 0.24], 0.04, 0.02);
OrPlots(3,:)   = getPlotHandles(1,2,[0.57 0.05 0.4 0.10], 0.04, 0.03);
PlotHs{1} = SizePlots;
PlotHs{2} = SfPlots;
PlotHs{3} = OrPlots;
PlotHs{4} = CnPlots;

for sub = 1:2
    subject = subjectNames{sub};
    [expDates, protocolNames, protocolTypes, electrodes, elecType] = getProtocolInfoGratings(subject);
    
    for p = 1:length(showForParams)
       par = showForParams{p};
       
       if strcmp(par, 'Size')   % for individual elecs
           clear BLpwr STpwr
           for el = 1:length(electrodes)
                elecis  = electrodes(el);
                expDate = expDates{2}{el}; 
                protocolName = protocolNames{2}{el};
                PSp_file = fullfile(dataDir,'derivatives','Spectra',['PowerSp_' subject expDate protocolName '_st_' num2str(1000*tST(1)) '_' num2str(1000*tST(2))... 
                    '_bl_' num2str(1000*tBL(1)) '_' num2str(1000*tBL(2)) '_method' num2str(methodFlag) '_notch' num2str(useNotchData) '.mat']);
                D = load(PSp_file);
                sizes   = D.parameterCombinations.sValsUnique; 
                sizesind= find(sizes>0.05);  % for some conditions
                paramX = sizes(sizesind);
                paramN = 'Size';
                if ~isempty(useoriforsze{sub})
                    useorind = find(D.parameterCombinations.oValsUnique==useoriforsze{sub});
                else useorind = [1:length(D.parameterCombinations.oValsUnique)];
                end
                if elecis>81
                    ind = D.EcogElectrodes == elecis;
                    BLpwr(el,:,:) = squeeze((permute(repmat(D.blPwrEcog(ind,:),[1 1 length(sizesind)]),[1 3 2])));
                    temp = permute(D.stPwrEcog(ind,sizesind,useorind,:),[2 3 4 1]);
                else
                    ind = D.LFPElectrodes == elecis;
                    BLpwr(el,:,:) = squeeze((permute(repmat(D.blPwrLfp(ind,:),[1 1 length(sizesind)]),[1 3 2])));
                    temp = permute(D.stPwrLfp(ind,sizesind,useorind,:),[2 3 4 1]);
                end
                % pool across ori conditions?
                for i = 1:size(temp,1)
                    temp2= [];
                    for j = 1:size(temp,2)
                        temp2 = cat(2,temp2,squeeze(temp(i,j,:)));
                    end
                    STpwr(el,i,:) = mean(temp2,2);
                end
           end
           lfpelecs = sum(electrodes<=81);
           ecogelecs= sum(electrodes>81);
           % lfp power % ecog power
           pwrDLfp = squeeze(nanmean(log10(STpwr((electrodes<=81),:,:)) - log10(BLpwr((electrodes<=81),:,:))));   % avg over electrodes
           pwrDEcog= squeeze(nanmean(log10(STpwr((electrodes>81),:,:))- log10(BLpwr((electrodes>81),:,:))));   % avg over electrodes
           % band power
           f_inds    = D.freqVals>=fBand(1) & D.freqVals<=fBand(2);
           stBandPwr = log10(squeeze(sum(STpwr(electrodes<=81,:,f_inds),3)));
           blBandPwr = log10(squeeze(sum(BLpwr(electrodes<=81,:,f_inds),3)));
           delPwrLfp = squeeze(nanmean((stBandPwr-blBandPwr),1));
                    
           stBandPwr = log10(squeeze(sum(STpwr(electrodes>81,:,f_inds),3)));
           blBandPwr = log10(squeeze(sum(BLpwr(electrodes>81,:,f_inds),3)));
           delPwrEcog= squeeze(nanmean((stBandPwr-blBandPwr),1));
           
       elseif strcmp(par, 'SF')
           clear BLpwr STpwr
           expDate = expDates{1}{1}; % same for all electrodes
           protocolName = protocolNames{1}{1};
           PSp_file = fullfile(dataDir,'derivatives','Spectra',['PowerSp_' subject expDate protocolName '_st_' num2str(1000*tST(1)) '_' num2str(1000*tST(2))...
                   '_bl_' num2str(1000*tBL(1)) '_' num2str(1000*tBL(2)) '_method' num2str(methodFlag) '_notch' num2str(useNotchData) '.mat']);
           D = load(PSp_file);
           paramN = 'SF';
           paramX   = D.parameterCombinations.fValsUnique;
           if ~isempty(useoriforsf{sub})
                useorind = find(D.parameterCombinations.oValsUnique==useoriforsf{sub});
           else useorind = [1:length(D.parameterCombinations.oValsUnique)];
           end 
           lfpelecs = size(D.stPwrLfp,1);
           ecogelecs = size(D.stPwrEcog,1);
           STpwr   = log10(squeeze(mean(D.stPwrLfp(:,:,useorind,:),3)));   % mean over ori conds
           BLpwr   = log10(permute(repmat(D.blPwrLfp,1,1,size(STpwr,2)),[1 3 2]));
           pwrDLfp = squeeze(mean(STpwr - BLpwr));   % avg over electrodes
           STpwr   = log10(squeeze(mean(D.stPwrEcog(:,:,useorind,:),3)));
           BLpwr   = log10(permute(repmat(D.blPwrEcog,1,1,size(STpwr,2)),[1 3 2]));
           pwrDEcog= squeeze(mean(STpwr - BLpwr));
           % band power
           f_inds    = D.freqVals>=fBand(1) & D.freqVals<=fBand(2);
           blBandPwr = log10(repmat(squeeze(sum(D.blPwrLfp(:,f_inds),2)),[1 length(paramX)]));
           stBandPwr = log10(squeeze(mean(sum(D.stPwrLfp(:,:,useorind,f_inds),4),3)));   % mean over oris
           delPwrLfp = squeeze(mean((stBandPwr-blBandPwr),1));
           blBandPwr = log10(repmat(squeeze(sum(D.blPwrEcog(:,f_inds),2)),[1 length(paramX)]));
           stBandPwr = log10(squeeze(mean(sum(D.stPwrEcog(:,:,useorind,f_inds),4),3)));
           delPwrEcog= squeeze(mean((stBandPwr-blBandPwr),1));
               
       elseif  strcmp(par, 'Ori') 
           clear BLpwr STpwr
           expDate = expDates{1}{1}; % same for all electrodes
           protocolName = protocolNames{1}{1};
           PSp_file = fullfile(dataDir,'derivatives','Spectra',['PowerSp_' subject expDate protocolName '_st_' num2str(1000*tST(1)) '_' num2str(1000*tST(2))...
                  '_bl_' num2str(1000*tBL(1)) '_' num2str(1000*tBL(2)) '_method' num2str(methodFlag) '_notch' num2str(useNotchData) '.mat']);
           D = load(PSp_file);
           paramN = 'Ori';
           paramX = D.parameterCombinations.oValsUnique;
           if ~isempty(usesfforori{sub})
                usesfind = find(D.parameterCombinations.fValsUnique==usesfforori{sub});  % sf condition to choose for ori tuning?
           else  usesfind = [1:length(D.parameterCombinations.fValsUnique)];
           end 
           
           lfpelecs = size(D.stPwrLfp,1);
           ecogelecs = size(D.stPwrEcog,1);
           STpwr   = log10(squeeze(mean(D.stPwrLfp(:,usesfind,:,:),2)));  % mean over sfs
           BLpwr   = log10(permute(repmat(D.blPwrLfp,1,1,size(STpwr,2)),[1 3 2]));
           pwrDLfp = squeeze(mean(STpwr - BLpwr));   % avg over electrodes
           STpwr   = log10(squeeze(mean(D.stPwrEcog(:,usesfind,:,:),2)));
           BLpwr   = log10(permute(repmat(D.blPwrEcog,1,1,size(STpwr,2)),[1 3 2]));
           pwrDEcog= squeeze(mean(STpwr - BLpwr));   % avg over electrodes
           % band power
           f_inds    = D.freqVals>=fBand(1) & D.freqVals<=fBand(2);
           blBandPwr = log10(repmat(squeeze(sum(D.blPwrLfp(:,f_inds),2)),[1 length(paramX)]));
           stBandPwr = log10(squeeze(mean(sum(D.stPwrLfp(:,usesfind,:,f_inds),4),2)));   % mean over sf s
           delPwrLfp = squeeze(mean((stBandPwr-blBandPwr),1));
           blBandPwr = log10(repmat(squeeze(sum(D.blPwrEcog(:,f_inds),2)),[1 length(paramX)]));
           stBandPwr = log10(squeeze(mean(sum(D.stPwrEcog(:,usesfind,:,f_inds),4),2)));
           delPwrEcog= squeeze(mean((stBandPwr-blBandPwr),1));
           
       elseif  strcmp(par, 'Con')
           clear BLpwr STpwr
           expDate = expDates{3}{1}; % same for all electrodes
           protocolName = protocolNames{3}{1};
           PSp_file = fullfile(dataDir,'derivatives','Spectra',['PowerSp_' subject expDate protocolName '_st_' num2str(1000*tST(1)) '_' num2str(1000*tST(2))...
                  '_bl_' num2str(1000*tBL(1)) '_' num2str(1000*tBL(2)) '_method' num2str(methodFlag) '_notch' num2str(useNotchData) '.mat']);
           D = load(PSp_file);
           paramN = 'Con';
           paramX = D.parameterCombinations.cValsUnique;
          if ~isempty(useoriforcon{sub})
              useorind = find(D.parameterCombinations.oValsUnique==useoriforcon{sub});
          else  useorind = [1:length(D.parameterCombinations.oValsUnique)];
          end 
           
           lfpelecs = size(D.stPwrLfp,1);
           ecogelecs = size(D.stPwrEcog,1);
           STpwr   = log10(squeeze(mean(D.stPwrLfp(:,:,useorind,:),3)));  % mean over oris
           BLpwr   = log10(permute(repmat(D.blPwrLfp,1,1,size(STpwr,2)),[1 3 2]));
           pwrDLfp = squeeze(mean(STpwr - BLpwr));   % avg over electrodes
           STpwr   = log10(squeeze(mean(D.stPwrEcog(:,:,useorind,:),3)));
           BLpwr   = log10(permute(repmat(D.blPwrEcog,1,1,size(STpwr,2)),[1 3 2]));
           pwrDEcog= squeeze(mean(STpwr - BLpwr));   % avg over electrodes
           % band power
           f_inds    = D.freqVals>=fBand(1) & D.freqVals<=fBand(2);
           blBandPwr = log10(repmat(squeeze(sum(D.blPwrLfp(:,f_inds),2)),[1 length(paramX)]));
           stBandPwr = log10(squeeze(mean(sum(D.stPwrLfp(:,:,useorind,f_inds),4),3)));
           delPwrLfp = squeeze(mean((stBandPwr-blBandPwr),1));
           blBandPwr = log10(repmat(squeeze(sum(D.blPwrEcog(:,f_inds),2)),[1 length(paramX)]));
           stBandPwr = log10(squeeze(mean(sum(D.stPwrEcog(:,:,useorind,f_inds),4),3)));
           delPwrEcog= squeeze(mean((stBandPwr-blBandPwr),1));
           
       else 
               disp('param not present'); return;
       end
       
       % plots:
       useplts = PlotHs{p};
       
       plot(useplts(1,sub),D.freqVals,zeros(size(D.freqVals)),'--','LineWidth',1,'color',[0.5 0.5 0.5]); hold(useplts(1,sub),'on');
       h = plot(useplts(1,sub),D.freqVals,10*(pwrDLfp),'LineWidth',1); 
       set(h,{'color'},num2cell(stimCols{p}(1:length(paramX),:),2)); 
       set(useplts(1,sub),'XTick',fTicks,'XTickLabel',[],'xlim',[0 150],'fontSize',fontSizeTicks,'fontweight','bold','tickLength',tickLength,'tickdir','out','box','off');
       
       plot(useplts(2,sub),D.freqVals,zeros(size(D.freqVals)),'--','LineWidth',1,'color',[0.5 0.5 0.5]); hold(useplts(2,sub),'on');
       h = plot(useplts(2,sub),D.freqVals,10*(pwrDEcog),'LineWidth',1); 
       set(h,{'color'},num2cell(stimCols{p}(1:length(paramX),:),2)); 
       set(useplts(2,sub),'XTick',fTicks,'XTickLabel',fTicks,'xlim',[0 150],'fontSize',fontSizeTicks,'fontweight','bold','tickLength',tickLength,'tickdir','out','box','off');
       
       hold(useplts(3,sub),'on');
       plot(useplts(3,sub),1:length(paramX),10*delPwrLfp, 'Color',lfpcol,'Marker','*');
       plot(useplts(3,sub),1:length(paramX),10*delPwrEcog,'Color',ecogcol,'Marker','*');
       set( useplts(3,sub),'xlim',[0.5,length(paramX)+0.5],'Xtick',[1:length(paramX)],'XtickLabel',round(paramX,1),'fontSize',fontSizeTicks,'fontweight','bold','tickLength',tickLength,'tickdir','out','box','off');
       
       xlabel(useplts(2,sub),'Frequency (Hz)','fontSize',fontSizeSmall,'fontWeight','bold');
       if sub==1
           text(-0.3,-0.7,'Change in power (dB)','Rotation',90,'Parent',useplts(1,sub),'Units','Normalized','fontSize',fontSizeSmall,'fontWeight','bold');  
           text(0.8,1.3,[paramN,' tuning'],'Parent',useplts(1,sub),'Units','Normalized','fontSize',fontSizeLarge,'fontWeight','bold');  
       end
       if sub==2
           delis = 0.1;
           for dd= 1:length(paramX)
               text(0.8,1-(dd-1)*delis,num2str(paramX(dd)),'Color',stimCols{p}(dd,:),'Parent',useplts(1,sub),'units','normalized','fontSize',fontSizeTicks,'fontWeight','bold')
           end
       end 
       text(0.80,0.1,['N=',num2str(lfpelecs)],'Parent',useplts(1,sub),'Color',lfpcol,'Units','Normalized','fontSize',fontSizeTicks,'fontWeight','bold')
       text(0.80,0.1,['N=',num2str(ecogelecs)],'Parent',useplts(2,sub),'Color',ecogcol,'Units','Normalized','fontSize',fontSizeTicks,'fontWeight','bold')
       title(useplts(1,sub),['M ',num2str(sub)],'fontSize',fontSizeSmall,'fontWeight','bold')
       
       xlabel(useplts(3,sub),paramN,'fontSize',fontSizeSmall,'fontWeight','bold');
       ylabel(useplts(3,1),'Band power (dB)','fontSize',fontSizeSmall,'fontWeight','bold');
       text(0.02,0.85,[num2str(fBand),' Hz'],'Parent',useplts(3,sub),'Units','Normalized','fontSize',fontSizeTicks,'fontWeight','bold')
       text(0.60,0.1,'LFP','Parent',useplts(3,sub),'Color',lfpcol,'Units','Normalized','fontSize',fontSizeTicks,'fontWeight','bold')
       text(0.80,0.1,'ECoG','Parent',useplts(3,sub),'Color',ecogcol,'Units','Normalized','fontSize',fontSizeTicks,'fontWeight','bold')
       
           
           
    end 
end 




%%%%
% for seeing the tuning for all conditions:

% figure;
% SizePlots = getPlotHandles(1,2,[0.07 0.54 0.4 0.20], 0.04, 0.03);
% SfPlots   = getPlotHandles(1,2,[0.57 0.54 0.4 0.20], 0.04, 0.03);
% CnPlots   = getPlotHandles(1,2,[0.07 0.05 0.4 0.20], 0.04, 0.03);
% OrPlots   = getPlotHandles(1,2,[0.57 0.05 0.4 0.20], 0.04, 0.03);
% PlotHs{1} = SizePlots;
% PlotHs{2} = SfPlots;
% PlotHs{3} = OrPlots;
% PlotHs{4} = CnPlots;
% 
% for sub = 1:2
%     subject = subjectNames{sub};
%     [expDates, protocolNames, protocolTypes, electrodes, elecType] = getProtocolInfoGratings(subject);
%     for p = 1:length(showForParams)
%        par = showForParams{p};
%        
%        if strcmp(par, 'Size')   % for individual elecs
%            clear BLpwr STpwr
%            for el = 1:length(electrodes)
%                 elecis  = electrodes(el);
%                 expDate = expDates{2}{el}; 
%                 protocolName = protocolNames{2}{el};
%                 PSp_file = fullfile(dataDir,'derivatives','Spectra',['PowerSp_' subject expDate protocolName '_st_' num2str(1000*tST(1)) '_' num2str(1000*tST(2))... 
%                     '_bl_' num2str(1000*tBL(1)) '_' num2str(1000*tBL(2)) '_method' num2str(methodFlag) '_notch' num2str(useNotchData) '.mat']);
%                 D = load(PSp_file);
%                 paramY  = D.parameterCombinations.oValsUnique; 
%                 sizes   = D.parameterCombinations.sValsUnique; 
%                 sizesind= find(sizes>0.05);  % for some conditions
%                 paramX = sizes(sizesind);
%                 paramN = 'Size';
%                 
%                 if elecis>81
%                     ind = D.EcogElectrodes == elecis;
%                     BLpwr(el,:,:,:) = squeeze((permute(repmat(D.blPwrEcog(ind,:),[1 1 length(sizesind) length(paramY)]),[1 3 4 2])));
%                     STpwr(el,:,:,:) = permute(D.stPwrEcog(ind,sizesind,:,:),[2 3 4 1]);
%                 else
%                     ind = D.LFPElectrodes == elecis;
%                     BLpwr(el,:,:,:) = squeeze((permute(repmat(D.blPwrLfp(ind,:),[1 1 length(sizesind) length(paramY)]),[1 3 4 2])));
%                     STpwr(el,:,:,:) = permute(D.stPwrLfp(ind,sizesind,:,:),[2 3 4 1]);
%                 end
%                 
%            end
%            lfpelecs = sum(electrodes<=81);
%            ecogelecs= sum(electrodes>81);
%            % band power
%            f_inds    = D.freqVals>=fBand(1) & D.freqVals<=fBand(2);
%            stBandPwr = log10(squeeze(sum(STpwr(electrodes<=81,:,:,f_inds),4)));
%            blBandPwr = log10(squeeze(sum(BLpwr(electrodes<=81,:,:,f_inds),4)));
%            delPwrLfp = squeeze(nanmean((stBandPwr-blBandPwr),1));  % size x ori
%                     
%            stBandPwr = log10(squeeze(sum(STpwr(electrodes>81,:,:,f_inds),4)));
%            blBandPwr = log10(squeeze(sum(BLpwr(electrodes>81,:,:,f_inds),4)));
%            delPwrEcog= squeeze(nanmean((stBandPwr-blBandPwr),1));
%            
%        elseif strcmp(par, 'SF')
%            clear BLpwr STpwr
%            expDate = expDates{1}{1}; % same for all electrodes
%            protocolName = protocolNames{1}{1};
%            PSp_file = fullfile(dataDir,'derivatives','Spectra',['PowerSp_' subject expDate protocolName '_st_' num2str(1000*tST(1)) '_' num2str(1000*tST(2))...
%                    '_bl_' num2str(1000*tBL(1)) '_' num2str(1000*tBL(2)) '_method' num2str(methodFlag) '_notch' num2str(useNotchData) '.mat']);
%            D = load(PSp_file);
%            paramN = 'SF';
%            paramX   = D.parameterCombinations.fValsUnique;
%            paramY = D.parameterCombinations.oValsUnique;
%            
%            lfpelecs = size(D.stPwrLfp,1);
%            ecogelecs= size(D.stPwrEcog,1);
%            % band power
%            f_inds    = D.freqVals>=fBand(1) & D.freqVals<=fBand(2);
%            blBandPwr = log10(repmat(squeeze(sum(D.blPwrLfp(:,f_inds),2)),[1 length(paramX) length(paramY)]));
%            stBandPwr = log10(squeeze(sum(D.stPwrLfp(:,:,:,f_inds),4)));   % mean over oris
%            delPwrLfp = squeeze(mean((stBandPwr-blBandPwr),1));
%            blBandPwr = log10(repmat(squeeze(sum(D.blPwrEcog(:,f_inds),2)),[1 length(paramX) length(paramY)]));
%            stBandPwr = log10(squeeze(sum(D.stPwrEcog(:,:,:,f_inds),4)));
%            delPwrEcog= squeeze(mean((stBandPwr-blBandPwr),1));    % sf x ori
%                
%        elseif  strcmp(par, 'Ori') 
%            clear BLpwr STpwr
%            expDate = expDates{1}{1}; % same for all electrodes
%            protocolName = protocolNames{1}{1};
%            PSp_file = fullfile(dataDir,'derivatives','Spectra',['PowerSp_' subject expDate protocolName '_st_' num2str(1000*tST(1)) '_' num2str(1000*tST(2))...
%                   '_bl_' num2str(1000*tBL(1)) '_' num2str(1000*tBL(2)) '_method' num2str(methodFlag) '_notch' num2str(useNotchData) '.mat']);
%            D = load(PSp_file);
%            paramN = 'Ori';
%            paramX = D.parameterCombinations.oValsUnique;
%            paramY = D.parameterCombinations.fValsUnique;
%            
%            lfpelecs = size(D.stPwrLfp,1);
%            ecogelecs = size(D.stPwrEcog,1);
%            % band power
%            f_inds    = D.freqVals>=fBand(1) & D.freqVals<=fBand(2);
%            blBandPwr = log10(repmat(squeeze(sum(D.blPwrLfp(:,f_inds),2)),[1 length(paramY) length(paramX)]));
%            stBandPwr = log10(squeeze(sum(D.stPwrLfp(:,:,:,f_inds),4)));   
%            delPwrLfp = squeeze(mean((stBandPwr-blBandPwr),1))';  % ori x sf
%            blBandPwr = log10(repmat(squeeze(sum(D.blPwrEcog(:,f_inds),2)),[1 length(paramY) length(paramX)]));
%            stBandPwr = log10(squeeze(sum(D.stPwrEcog(:,:,:,f_inds),4)));
%            delPwrEcog= squeeze(mean((stBandPwr-blBandPwr),1))';
%            
%        elseif  strcmp(par, 'Con')
%            clear BLpwr STpwr
%            expDate = expDates{3}{1}; % same for all electrodes
%            protocolName = protocolNames{3}{1};
%            PSp_file = fullfile(dataDir,'derivatives','Spectra',['PowerSp_' subject expDate protocolName '_st_' num2str(1000*tST(1)) '_' num2str(1000*tST(2))...
%                   '_bl_' num2str(1000*tBL(1)) '_' num2str(1000*tBL(2)) '_method' num2str(methodFlag) '_notch' num2str(useNotchData) '.mat']);
%            D = load(PSp_file);
%            paramN = 'Con';
%            paramX = D.parameterCombinations.cValsUnique;
%            paramY = D.parameterCombinations.oValsUnique;
%            
%            lfpelecs = size(D.stPwrLfp,1);
%            ecogelecs = size(D.stPwrEcog,1);
%            % band power
%            f_inds    = D.freqVals>=fBand(1) & D.freqVals<=fBand(2);
%            blBandPwr = log10(repmat(squeeze(sum(D.blPwrLfp(:,f_inds),2)),[1 length(paramX) length(paramY)]));
%            stBandPwr = log10(squeeze(sum(D.stPwrLfp(:,:,:,f_inds),4)));
%            delPwrLfp = squeeze(mean((stBandPwr-blBandPwr),1));  % con x ori
%            blBandPwr = log10(repmat(squeeze(sum(D.blPwrEcog(:,f_inds),2)),[1 length(paramX) length(paramY)]));
%            stBandPwr = log10(squeeze(sum(D.stPwrEcog(:,:,:,f_inds),4)));
%            delPwrEcog= squeeze(mean((stBandPwr-blBandPwr),1));
%            
%        else 
%                disp('param not present'); return;
%        end
%        
%        % plots:
%        useplts = PlotHs{p};
%        
%        hold(useplts(1,sub),'on');
%        h = plot(useplts(1,sub),1:length(paramX),10*delPwrLfp,'Marker','*');
%        set(h,{'color'},num2cell(stimCols{3}(1:length(paramY),:),2));
%        h = plot(useplts(1,sub),1:length(paramX),10*delPwrEcog,'Marker','o');
%        set(h,{'color'},num2cell(stimCols{3}(1:length(paramY),:),2));
%        set( useplts(1,sub),'xlim',[0.5,length(paramX)+0.5],'Xtick',[1:length(paramX)],'XtickLabel',round(paramX,1),'fontSize',fontSizeTicks,'fontweight','bold','tickLength',tickLength,'tickdir','out','box','off');
%        
%        title(useplts(1,sub),['M ',num2str(sub)],'fontSize',fontSizeSmall,'fontWeight','bold')
%        xlabel(useplts(1,sub),paramN,'fontSize',fontSizeSmall,'fontWeight','bold');
%        ylabel(useplts(1,1),'Band power (dB)','fontSize',fontSizeSmall,'fontWeight','bold');
%        text(0.02,0.85,[num2str(fBand),' Hz'],'Parent',useplts(1,sub),'Units','Normalized','fontSize',fontSizeTicks,'fontWeight','bold') 
%        
%        if sub==2
%            delis = 0.1;
%            for dd= 1:length(paramY)
%                text(0.8,1-(dd-1)*delis,num2str(paramY(dd)),'Color',stimCols{3}(dd,:),'Parent',useplts(1,sub),'units','normalized','fontSize',fontSizeTicks,'fontWeight','bold')
%            end
%        end 
%            
%     end 
% end

    
