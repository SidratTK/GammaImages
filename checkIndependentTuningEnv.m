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
            if strcmp(protocolTypes{pt},'SizeOri') || elecis>81 
                indel = find(electrodes_sz==elecis);
                pName =  protocolNames{psz}{indel}; eDate =  expDates{psz}{indel};
            else
                pName =  ''; eDate =  '';
            end
            namedetail = ['dGamma_',namedetail1,num2str(elecis),'_',protocolTypes{1},'_',protocolTypes{2},eDate,pName,'_',protocolTypes{3},namedetail2,'.mat'];
            fname = fullfile(dataDir,'derivatives','dGamma',namedetail);
            dG    = load(fname);
            delGuse = squeeze(dG.delGamma{pt});
            [r,ptemp,g,sepind] = checkIndependentTuning(delGuse);
            if strcmp(subject,'alpaH') && elecis==41     % display for typical elec
                varH = dG.parameterCombinations{pt}.oValsUnique;
                if strcmp(protocolTypes{pt},'SizeOri'),
                    varV = dG.parameterCombinations{pt}.sValsUnique;
                elseif strcmp(protocolTypes{pt},'SFOri')
                    varV = dG.parameterCombinations{pt}.fValsUnique;
                elseif strcmp(protocolTypes{pt},'ConOri')
                    varV = dG.parameterCombinations{pt}.cValsUnique;
                end
                [r,ptemp,g,sepind] = checkIndependentTuning(delGuse,1,varH,varV);
            end
            rsq{sub,pt}(e) = r^2;
            p{sub,pt}(e)   = ptemp; 
            gof{sub,pt}(e) = g;
            sep{sub,pt}(e) = sepind;
        end 
    end
end
printf('\n');
printf('Squared correlation R2 represents how much variance is shared between 2 distributions'); 
printf('Or how much variance of A(observed-data) can be explained by B(separable-tuning-estimate)')
printf('GOF is a goodness of fit gives the error metric between the two distributions like MSE');
printf('Separability index is calculated from SVD decomposition and tells us how well the data can be factored into independent components');

palph = 0.05;
for sub = 1:2
    printf(subjectNames{sub});
    for pt=1:3
        r2 = round(mean(rsq{sub,pt}),2);    sr = round(std(rsq{sub,pt}),2);
        n  = length(rsq{sub,pt});
        s  = sum(p{sub,pt} <= palph);
        g  = round(mean(gof{sub,pt}),2);    sg = round(std(gof{sub,pt}),2);
        sepind = round(mean(sep{sub,pt}),2);ss = round(std(sep{sub,pt}),2);
        
        disp([' ', protocolTypes{pt},':',', n = ',num2str(n)]);
        printf(['\t r2 = ',num2str(r2),'+-',num2str(sr), 'SD. Significant p-val for ',num2str(s),' of ',num2str(n)]);
        printf(['\t GOF = ',num2str(g),'+-',num2str(sg), 'SD.\n\t Separability index = ',num2str(sepind),'+-',num2str(ss),'SD']);
    end
end


end