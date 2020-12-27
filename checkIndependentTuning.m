% this code is to check independence
% Responses to 2 variables are input 
% like SFOri or SizeOri data. responses can be LFP power in gamma band.
% computes marginal responses by averaging across rows or cols
% computes estimated joint matrix if 2 distributions were independent- 
% product of marginals
% see also Mazer et al., 2002 PNAS who follow a similar marginal analysis
% the observed joint tuning (data) and estimated separable tuning (prod of
% marginals) are compared via the Pearson correlation coefficient, the square of
% which represents the fraction of variance in the data shared by the estimate.
% Note that the product of the marginal tuning may not be the same order of
% magnitude as the observed data, but the correl coeff is not affected by
% intercept. 
% also calculate separability index as in Mazer et al.
% see Pena & Konishi compare additive (marginals) & multiplicative (factors) models
% use Model fitting to get the separable estimates. The R2 in linear model is
% v.v.close to vals obtained by data & importantly gives F statistic.
% REFS: Mazer et al., 2002, PNAS, Pena & Konishi, 2001, Science,

% STK 091220
% 101220: normalise data2D so it behaves like a prob distribn. Make sure only +ve values are input
% 171220: added option of plotting marginals 
% 221220: included SVD. and model estimates
% 231220: add the coef test to get f stat of models. change the output format

% make dummy data example
% data2DIn = rand(4,3);        % 2 dimensional random data OR
% parV = [1:8]; mparV = 4; sparV = 2; margV1 = pdf('Normal',parV,mparV,sparV); margV1 = margV1(:); 
% parH = [1:4]; mparH = 2; sparH = 1; margH1 = pdf('Normal',parH,mparH,sparH);
% data2DIn = margV1 * margH1;
% checkIndependentTuning(data2DIn,1)

function [SVDPrdOut,sepindex,MargPrdOut,MargSumOut,MargBthOut,chi2MargPrdOut,mrgsvdcorr] = checkIndependentTuning(data2DIn,showfig,labH,labV)
if ~exist('showfig','var'), showfig=0;   end

varH   = 1:size(data2DIn,2);
varV   = 1:size(data2DIn,1);
if ~exist('labH','var'), labH = varH; end
if ~exist('labV','var'), labV = varV; end

% get marginal sums of row & col
sumVa = repmat(sum(data2DIn,2),[1,size(data2DIn,2)]);
sumHa = repmat(sum(data2DIn,1),[size(data2DIn,1),1]);

% data values tells how likely a chosen input combination is to generate a strong response. 
% if want the data to behave like probabilities, sum should be 1, so normalise
% the stimulus combs shown independently, => 2D data is independently sampled. 
% the 3 conditions of being a pdf are met: +ve, sum=1, independent samples 
sumis  = sum(sum(data2DIn));
data2D = data2DIn/sumis;  % subsequent linear corrs hold even if this normalisation not done.

% 1. get marginals (normalised)
margV = sum(data2D,2);    % marginal of Vertical param
margH = sum(data2D,1);    % marginal of Horizontal param

% 2. get SVD factors,
[U,S,V] = svd(data2DIn,'econ');  
S = diag(S); factH = V(:,1); factV = U(:,1);    % val & first singular vectors
sepindex = S(1)^2/(sum(S.^2));

% 3. correlation between the factors
[rh,ph] = corr(margH', factH);
[rv,pv] = corr(margV , factV);

% 4A. separable product of marginals
data2DIndp = margV*margH*sumis;     % column vector x row vector
[r,pr] = corr([data2DIndp(:), data2DIn(:)]);  % Corr coeff statistic as in Mazer et al 2002, for R2
rdta_margprd = r(2); pdta_margprd = pr(2);

% 4B. separable sum of marginals
data2DInds  = sumis*(repmat(margV,[1,size(data2D,2)]) + repmat(margH,[size(data2D,1),1]));  % col val + row val
[r,pr]      = corr([data2DInds(:), data2DIn(:)]);
rdta_margsum= r(2); pdta_margsum = pr(2);

% 4C. separable product of svd factors 
data2Dsvd = factV * S(1,1) * factH';      
[r,pr]    = corr([data2Dsvd(:), data2DIn(:)]);  % Correlation coeff for getting R2
rdta_svdprd = r(2); pdta_svdprd = pr(2);

% 5. calc the chi2 between the dataIn and product of marginals, 
wid = 'stats:chi2gof:LowCounts';   warning('off',wid);
bincntrs = 1:length(varH)*length(varV);   %  get bin numbers
[hchi,pchi,st] = chi2gof(bincntrs,'Ctrs',bincntrs,'Frequency',data2DIn(:),'Expected',data2DIndp(:));   
warning('on',wid);

% 6A. use model fits for the above marginal model & svd model
% use data2DIn vs data2D or sumH vs margH doesnt change the signficance but only model coeff magnitudes
y    = data2DIn(:);
sumH = sumHa(:);   
sumV = sumVa(:);
d    = dataset(y, sumH, sumV);
LMsum_marg = LinearModel.fit(d,'y ~ sumH + sumV');
LMprd_marg = LinearModel.fit(d,'y ~ sumH : sumV');
LMbth_marg = LinearModel.fit(d,'y ~ sumH * sumV');  % just check both interactions
fH = sqrt(S(1))* repmat(factH(:)',[size(data2DIn,1),1]); fH = fH(:);
fV = sqrt(S(1))* repmat(factV(:) ,[1,size(data2DIn,2)]); fV = fV(:);
d2 = dataset(y, fH, fV);
LMprd_svd = LinearModel.fit(d2,'y ~ fH : fV');

% 6B. pick up some features from these models here:
r2svdprd = LMprd_svd.Rsquared.Ordinary;   % R2 reps ratio of variance explained by the model
r2margprd= LMprd_marg.Rsquared.Ordinary;
r2margsum= LMsum_marg.Rsquared.Ordinary;
r2margbth= LMbth_marg.Rsquared.Ordinary;
[p_svdprd, f_svdprd, d_svdprd]  = coefTest(LMprd_svd); % p<palpha - model significantly diff than intercept-only model
[p_margprd,f_margprd,d_margprd] = coefTest(LMprd_marg); % f = f value from the f test
[p_margsum,f_margsum,d_margsum] = coefTest(LMsum_marg);
[p_margbth,f_margbth,d_margbth] = coefTest(LMbth_marg);

% 7. Update the estimates using the models
data2DIndp = predict(LMprd_marg,[sumH,sumV]); data2DIndp = reshape(data2DIndp,size(data2DIn));
data2DInds = predict(LMsum_marg,[sumH,sumV]); data2DInds = reshape(data2DInds,size(data2DIn));
data2Dsvd  = predict(LMprd_svd ,[fH,fV]);     data2Dsvd  = reshape(data2Dsvd ,size(data2DIn));

% 8. display
if showfig
    figH = figure;
    figH.Units = 'centimeters';
    figH.PaperType  = 'a4';
    figH.PaperUnits = 'centimeters';
    figH.PaperSize  = [17.5 15];
    figH.PaperOrientation = 'Portrait';
    figH.PaperPosition = [0 0 figH.PaperSize];
    figH.Color = [1 1 1];
    figH.Position = [0 0 figH.PaperSize];
    figH.PaperUnits = 'normalized';
    
    plotMh = getPlotHandles(1,1,[0.22 0.80 0.30 0.15], 0.02,0.02);
    plots1 = getPlotHandles(1,1,[0.61 0.80 0.30 0.15], 0.02,0.02);
    plotMv = getPlotHandles(1,1,[0.01 0.45 0.10 0.25], 0.02,0.02);
    plot2D = getPlotHandles(1,1,[0.22 0.45 0.30 0.25], 0.02,0.02);
    plot2Dsum = getPlotHandles(1,1,[0.22 0.10 0.30 0.25], 0.02,0.02);
    plot2Dprd = getPlotHandles(1,1,[0.61 0.10 0.30 0.25], 0.02,0.02);    
    plot2Dsvd = getPlotHandles(1,1,[0.61 0.45 0.30 0.25], 0.02,0.02);
    
    cols = gray(10);
    hold(plotMh,'on');  hold(plotMv,'on'); camroll(plotMv,90);
%     t = plot(plotMh,varH,data2DIn./sumVa,'LineWidth',1,'Linestyle',':');
%     set(t,{'color'},num2cell(cols(1:length(varV),:),2));
%     t = plot(plotMv,varV,data2DIn./sumHa,'LineWidth',1,'Linestyle',':');
%     set(t,{'color'},num2cell(cols(1:length(varH),:),2));
    
    bar(plots1,[1:length(S)],100*(sqrt(S.^2)./sum(sqrt(S.^2))));   % magnitude of sing vectors
    
    axes(plot2D); imagesc(flip(data2DIn,1));
    title(plot2D,'2D Tuning Data Y');
    
    axes(plot2Dprd); imagesc(flip(data2DIndp,1));
    title(plot2Dprd,'Y ~ 1 + F * G');
    
    axes(plot2Dsum); imagesc(flip((data2DInds),1));
    title(plot2Dsum,'Y ~ 1 + F + G');
    
    axes(plot2Dsvd); imagesc(flip(data2Dsvd,1));
    title(plot2Dsvd,'Y ~ 1 + U(1)*s1*V(1)');
    
    plot(plotMh,varH,margH,'LineWidth',2,'Color','b'); 
    plot(plotMv,varV,margV,'LineWidth',2,'Color','b'); 
    
    plot(plotMh,varH,factH/sum(factH),'LineWidth',2,'Color','r','Linestyle','--');
    plot(plotMv,varV,factV/sum(factV),'LineWidth',2,'Color','r','Linestyle','--'); 
    
    set(plot2D,'XTick',varH,'XTickLabel',labH,'YTick',varV,'YTickLabel',flip(labV),'XTickLabelRotation',0,'tickdir','out','box','off','fontWeight','bold');
    set(plot2Dprd,'XTick',varH,'XTickLabel',labH,'YTick',varV,'YTickLabel',flip(labV),'XTickLabelRotation',0,'tickdir','out','box','off','fontWeight','bold')
    set(plot2Dsum,'XTick',varH,'XTickLabel',labH,'YTick',varV,'YTickLabel',flip(labV),'XTickLabelRotation',0,'tickdir','out','box','off','fontWeight','bold')
    set(plot2Dsvd,'XTick',varH,'XTickLabel',labH,'YTick',varV,'YTickLabel',flip(labV),'XTickLabelRotation',0,'tickdir','out','box','off','fontWeight','bold')
    
    set(plotMh,'XTick',varH,'XTickLabel',labH,'tickdir','out','box','off','fontWeight','bold');
    set(plotMv,'XTick',varV,'XTickLabel',labV,'tickdir','out','box','off','fontWeight','bold');  
    
    set(plots1, 'box','off');
    
    text(1.01,0.9, ['R2=',num2str(round(r2margprd,2))],'Parent',plot2Dprd,'Units','normalized','FontSize',10);
    text(1.01,0.9, ['R2=',num2str(round(r2svdprd,2)) ],'Parent',plot2Dsvd,'Units','normalized','FontSize',10);
    text(-0.4,0.9, ['R2=',num2str(round(r2margsum,2))],'Parent',plot2Dsum,'Units','normalized','FontSize',10);
%     text(1.01,0.8, 'Ftest','Parent',plot2Dprd,'Units','normalized','FontSize',10);
%     text(1.01,0.8, 'Ftest','Parent',plot2Dsvd,'Units','normalized','FontSize',10);
%     text(-0.1,0.8, 'Ftest','Parent',plot2Dsum,'Units','normalized','FontSize',10);
%     text(1.01,0.7, ['p =',num2str(round(p_margprd,2))],'Parent',plot2Dprd,'Units','normalized','FontSize',10);
%     text(1.01,0.7, ['p =',num2str(round(p_svdprd,2))], 'Parent',plot2Dsvd,'Units','normalized','FontSize',10);
%     text(-0.1,0.7, ['p =',num2str(round(p_margsum,2))],'Parent',plot2Dsum,'Units','normalized','FontSize',10);

    printf('\n');
    disp(['Variance explained by product of marginals is ',num2str(100*(rdta_margprd.^2)),'% ']);
    disp(['Variance explained by model fit of marginals product is ',num2str(100*r2margprd),'% .',...
        ' model F-val= ',num2str(round(f_margprd,2)),' p= ',num2str(round(p_margprd,3))]);
    disp(['Chi2 test for difference bw observed & independent model data is h = ',num2str(hchi)]);
    disp(['Variance explained by sum of marginals is ',num2str(100*(rdta_margsum.^2)),'% ']);
    disp(['Variance explained by model fit of marginals sum is ',num2str(100*r2margsum),'% .',...
        ' model F-val= ',num2str(round(f_margsum,2)),' p= ',num2str(round(p_margsum,3))]);
    disp(['Separability index from SVD factors is ',num2str(sepindex)]);
    disp(['Variance explained by product of 1st SVD vectors is ',num2str(100*(rdta_svdprd.^2)),'% ']);
    disp(['Variance explained by model fit of svd product is ',num2str(100*r2svdprd),'% .',...
        ' model F-val= ',num2str(round(f_svdprd,2)),' p= ',num2str(round(p_svdprd,3))]);
end

% 9. compile outputs
chi2MargPrdOut = [hchi,pchi,st.chi2stat];
MargPrdOut.Mdl = LMprd_marg; MargPrdOut.r2=r2margprd; MargPrdOut.fstat=[f_margprd,p_margprd,d_margprd];
MargSumOut.Mdl = LMsum_marg; MargSumOut.r2=r2margsum; MargSumOut.fstat=[f_margsum,p_margsum,d_margsum];
MargBthOut.Mdl = LMbth_marg; MargBthOut.r2=r2margbth; MargBthOut.fstat=[f_margbth,p_margbth,d_margbth];
SVDPrdOut.Mdl  = LMprd_svd;  SVDPrdOut.r2 =r2svdprd;  SVDPrdOut.fstat =[f_svdprd, p_svdprd, d_svdprd ];
mrgsvdcorr     = [rh rv ph pv];

end
