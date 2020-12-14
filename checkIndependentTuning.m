% this code is to check independence
% Responses to 2 variables are input 
% like SFOri or SizeOri data. responses can be LFP power in gamma band.
% computes marginal responses by averaging across rows or cols
% computes estimated joint matrix if 2 distributions were independent- 
% product of marginals
% see also Mazer et al., 2002 PNAS who follow a similar marginal analysis
% the observed joint tuning (data) and estimated separable tuning (prod of
% marginals) are compared via the Pearson correlation coefficient, the of
% which represents the fraction of variance in the data shared by the
% estimate.
% Note that the product of the marginal tuning may not be the same order of
% magnitude as the observed data, but the correl coeff is not affected by
% scaling or shifting. So it is the relative trend in the 2d data that is
% being compared.
% STK 091220
% 101220: updates to normalise data2D so it behaves like a prob distribn.
% Make sure only +ve values are input


% make dummy data example
% data2D = rand(3,3);        % 2 dimensional data
% parV = [0:pi/4:pi]; mparV = pi/2; sparV = pi/8; margV1 = pdf('Normal',parV,mparV,sparV); margV1 = margV1(:); 
% parH = [1:4];     mparH = 2; sparH = 1;    margH1 = pdf('Normal',parH,mparH,sparH);
% data2DIn = margV1 * margH1;

function [rout,pout,gof,sepindex] = checkIndependentTuning(data2DIn,showfig)
if ~exist('showfig','var'), showfig=0;   end

% data values tells how likely a chosen input combination is to generate a strong response. 
% want the data to behave like probabilities, so sum should be 1
% as the stimulus combinations have been shown independently, it means each
% element of 2D data is independently sampled. 
% the 3 conditions of being a pdf are met: +ve, sum=1, independent samples

data2D = data2DIn/sum(sum(data2DIn));

% get marginals
margV = sum(data2D,2);    % marginal of Vertical param
margH = sum(data2D,1);    % marginal of Horizontal param

% separable product
data2DInd = margV*margH;   % column vector x row vector

% try with GOF, calculates error values
errtype = 'MSE';
gof = goodnessOfFit(data2DInd(:),data2D(:),errtype);   % 

% r2 statistic as in Mazer et al 2002 Pnas
[r,p] = corr([data2DInd(:), data2D(:)]);
rout = r(2); pout=p(2);

% the SVD route  % also in Mazer et al, similar to Depireux et al. 2001 JNP
[U,S,V] = svd(data2DIn,'econ');
singvec=diag(S);
sepindex = singvec(1)^2/(sum(singvec.^2));

% display
if showfig
    figure;
    subplot(211);
    imagesc(data2D);
    title('2D tuning');
    subplot(212);
    imagesc(data2DInd);
    title('2D separable tuning');
    
    printf('\n');
    disp(['Variance explained by independent model is ',num2str(100*(r(2).^2)),'% ']);
    disp(['Separability index from SVD factors is ',num2str(sepindex)]);
    disp([errtype,' GoodnessOfFit of independent model is ',num2str(gof)]);
end

end
