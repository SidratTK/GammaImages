% envelope function for checking independence of features
% this function takes all electrodes of the protocol type chosen and for
% each of them gets the independence between the variable checked 

% Output is value of separability of variables for all protocol types
% for both monkeys

% added option of plotting for 1 elec-171220

function [] = checkIndependentTuningEnv()

% these paths may be used
rootPath = gammaModelPath_st();
dataDir  = fullfile(rootPath, 'Data');

% power related
useNotchData = 0;
PwrmethodFlag= 0; % 0= MT & dB change in gamma, 1= Gaussian gamma fit like Hermes etal.
tBL = [-0.25 0];
tST = [0.25 0.5];
numTapers = 3;
fBandGamma= [35 70];

subjectNames = {'alpaH','kesariH'};
[~,~,LFPElectrodeList,EcogElectrodeList,~] = getRFdetails(subjectNames,dataDir);

for sub = 1:2
    subject = subjectNames{sub};
    [expDates, protocolNames, protocolTypes, electrodes_sz, electrodeTypes] = getProtocolInfoGratings(subject);
    psz = strcmp(protocolTypes,'SizeOri');
    namedetail1 = [num2str(fBandGamma(1)),'_',num2str(fBandGamma(2)),'_tST_',num2str(tST(1)*1000),'_',num2str(tST(2)*1000),'_Tapers',num2str(numTapers),'_',subject,'_elec'];
    namedetail2 = ['_method',num2str(PwrmethodFlag),'_notch',num2str(useNotchData)];
    
    for pt = 1:length(protocolTypes)
        
        if strcmp(protocolTypes{pt},'SizeOri')       
            electrodes = electrodes_sz;      % use only size protocol elecs
        else                                 % else use all good elecs.
            electrodes = cat(1, EcogElectrodeList{sub}, LFPElectrodeList{sub});   
        end 
        
        for e = 1:length(electrodes)
            elecis = electrodes(e);
            pName =  ''; eDate =  '';
            if strcmp(protocolTypes{pt},'SizeOri') || elecis>81 
                indel = find(electrodes_sz==elecis);
                pName = protocolNames{psz}{indel}; eDate =  expDates{psz}{indel};
            end
            namedetail = ['dGamma_',namedetail1,num2str(elecis),'_',protocolTypes{1},'_',protocolTypes{2},eDate,pName,'_',protocolTypes{3},namedetail2,'.mat'];
            fname = fullfile(dataDir,'derivatives','dGamma',namedetail);
            dG    = load(fname);
            delGuse = squeeze(dG.delGamma{pt});
            
            varH = dG.parameterCombinations{pt}.oValsUnique;
            if strcmp(protocolTypes{pt},'SizeOri'),
                varV = dG.parameterCombinations{pt}.sValsUnique;
            elseif  strcmp(protocolTypes{pt},'SFOri')
                varV = dG.parameterCombinations{pt}.fValsUnique;
            elseif  strcmp(protocolTypes{pt},'ConOri')
                varV = dG.parameterCombinations{pt}.cValsUnique;
            end 
            
            if strcmp(subject,'alpaH') && elecis==41     % display for typical elec
                showflg = 1;
            else
                showflg = 0;
            end
            
            [SVDPrd{sub,pt}{e},sepind,MrgPrd{sub,pt}{e},MrgSum{sub,pt}{e},MrgBth{sub,pt}{e},chi2MrgPrd,mrgsvdcorr] ...
                = checkIndependentTuning(delGuse,showflg,varH,varV);
            
            sep{sub,pt}(e) = sepind;
            hchi{sub,pt}(e)= chi2MrgPrd(1);
            rmguv{sub,pt}(e,:)= mrgsvdcorr;
            
        end 
    end
end
printf('\n');
printf('Squared correlation R2 represents how much variance is shared between 2 distributions'); 
printf('Or how much variance of A(observed-data) can be explained by B(separable-tuning-estimate)')
printf('Separability index is calculated from SVD decomposition and tells us how well the data can be factored into independent components');
printf('h-value of chi2test tests whether the observed values could have been drawn from the expected distribution');

palph = 0.05;
for sub = 1:2
    printf(subjectNames{sub});
    for pt=1:3
        
        sepind = round(mean(sep{sub,pt}),2);ss = round(std(sep{sub,pt}),2);
        his    = sum(hchi{sub,pt} == 0);
        
        % FOR NOW DISPLAYING ALL :: 
        r = cell2mat(cellfun(@(x) x.r2, SVDPrd{sub,pt},'un',0));
        r2svd = round(mean(r),2); sr2svd = round(std(r),2);
        n     = length(r);
        pfsvd = cell2mat(cellfun(@(x) x.fstat(2), SVDPrd{sub,pt},'un',0));
        pfsvd = sum(pfsvd <= palph);
        
        r1 = cell2mat(cellfun(@(x) x.r2, MrgPrd{sub,pt},'un',0));
        r2prd = round(mean(r1),2); sr2prd = round(std(r1),2);
        pfprd = cell2mat(cellfun(@(x) x.fstat(2), MrgPrd{sub,pt},'un',0));
        pfprd = sum(pfprd <= palph);
        
        r2 = cell2mat(cellfun(@(x) x.r2, MrgSum{sub,pt},'un',0));
        r2sum = round(mean(r2),2); sr2sum = round(std(r2),2);
        pfsum = cell2mat(cellfun(@(x) x.fstat(2), MrgSum{sub,pt},'un',0));
        pfsum = sum(pfsum <= palph);
        
        r3 = cell2mat(cellfun(@(x) x.r2, MrgBth{sub,pt},'un',0));
        r2bth = round(mean(r3),2); sr2bth = round(std(r3),2);
        pfbth = cell2mat(cellfun(@(x) x.fstat(2), MrgBth{sub,pt},'un',0));
        pfbth = sum(pfbth <= palph);
        
        disp([' ', protocolTypes{pt},':',', n = ',num2str(n)]);
        printf(['\t r2 = ',num2str(r2svd),'+-',num2str(sr2svd), 'SD. Significant F-stat p-val for ',num2str(pfsvd),' of ',num2str(n)]);
        printf(['\t Separability index (svd) = ',num2str(sepind),'+-',num2str(ss),'SD']);
        
        printf(['\t r2 = ',num2str(r2prd),'+-',num2str(sr2prd), 'SD. Significant F-stat p-val for ',num2str(pfprd),' of ',num2str(n)]);
        printf(['\t h = 0 (fail to reject same distribution) in , ', num2str(his),' of ',num2str(n)]);
    
        printf(['\t r2 = ',num2str(r2sum),'+-',num2str(sr2sum), 'SD. Significant F-stat p-val for ',num2str(pfsum),' of ',num2str(n)]);

        printf(['\t r2 = ',num2str(r2bth),'+-',num2str(sr2bth), 'SD. Significant F-stat p-val for ',num2str(pfbth),' of ',num2str(n)]);

        [h,p] = ttest(r2,r1); % t-test for whether R2 of MrgPrd is significantly diff from MrgSum
        printf(['\t paired ttest on R2 of Marg Sum & Prd. H = ',num2str(h),'P = ',num2str(round(p,3))]);
        
        inds = (rmguv{sub,pt}(:,3)<palph) & (rmguv{sub,pt}(:,4)<palph);
        mr = mean(abs([rmguv{sub,pt}(inds,1);rmguv{sub,pt}(inds,2)]));
        sr = std(abs([rmguv{sub,pt}(inds,1);rmguv{sub,pt}(inds,2)]));
        printf(['\t Corr bw SVD factors and marginals was significant for ',num2str(sum(inds)),' of ',num2str(n)]);
        printf(['\t Avg mag of Corr coef =', num2str(round(mr,2)),' +-',num2str(round(sr,2)),' SD']);
        
    end
end


end