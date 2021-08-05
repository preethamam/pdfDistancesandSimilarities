function [distance, similarity, corr_coeffient] = pdfsDistanceSimilarity(P, Q, method, p)
%%***********************************************************************%
%*                 PDF Distances and similarities                       *%
%* Probability density functions distances and similarities for various *%
%* family of distances functions.                                       *%
%* Code authors: Preetham Manjunatha and Liu Wu                         *%
%* Date: 07/05/2021                                                     *%
%************************************************************************%
%
%************************************************************************%
% Citation: Cha, Sung-Hyuk. "Comprehensive survey on distance/similarity 
% measures between probability density functions." City 1, no. 2 (2007): 1.
%
% Usage: [distance, similarity] = pdfsDistanceSimilarity(P, Q, method, p)
% 
% Inputs:  P          - Probability Density Function 1
%          Q          - Probability Density Function 2
%          method     - Distance method (e.g. "Euclidean"). 
%                       Available:
%                       'Euclidean', "CityBlock", "Minkowski", "Chebyshev", 
%                       "Sorensen", "Gower", "Soergel", "Kulczynski", "Canberra", 
%                       "Lorentzian", 'Intersection', 'WaveHedges','Czekanowski', 
%                       'Motyka', 'Kulczynski', 'Ruzicka', 'Tanimoto', 
%                       "InnerProduct", "HarmonicMean", "Cosine", "Kumar-Hassebrook", 
%                       "Jaccard", "Dice", 'Fidelity', 'Bhattacharya', 'Hellinger',
%                       'Matusita', 'Squared-chord', 'SquaredEuclidean', 'Pearson', 
%                       'Neyman', 'Squared', 'ProbabilisticSymmetrix', 'Divergence', 
%                       'Clark', 'AdditiveSymmetric', 'Kullback-Leibler', 'Jeffreys', 
%                       'KDivergence', 'Topsoe', 'Jensen-Shannon', 'JensenDifference', 
%                       'Taneja', 'Kumar-Johnson', 'Average'
%          p          - Power term for Minkowski distance method only.
%                       Default: 2
%
% Outputs: distance         - Distance between pdf 1 and pdf 2
%          similarity       - Similarity between pdf 1 and pdf 2
%          corr_coeffient   - Correlation coefficient between pdf 1 and pdf 2
% 
%  It is highly recommended to read the cited reference. Note: some of the
%  distances and similarities are NaN due to the division by zero and the 
%  log of zero. These cases deserve special attention.
%  See also: histcounts, ksdensity

% Inputs check
if nargin < 3
    method = 'Euclidean';
end

if strcmp(method, "Minkowski") && nargin < 4
    p = 2;
end

% Distances and similarities
switch method
    %----------------------------------------------------------------------
    % Lp Minkowski family
    case "Euclidean"
        distance = sqrt(sum(abs(P-Q).^2));
        similarity = 1 - distance;

    case "CityBlock"
        distance = sum(abs(P-Q));
        similarity = nan;

    case "Minkowski"
        distance = realpow(sum(abs(P-Q).^p), 1/p);
        similarity = 1 - distance;

    case "Chebyshev"
        distance = max(abs(P-Q));
        similarity = 1 - distance;
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % L1 family
    case "Sorensen"
        distance = sum(abs(P-Q))/sum(P+Q);
        similarity = 1 - distance;

    case "Gower"
        distance = sum(abs(P-Q))/length(P);
        similarity = 1 - distance;

    case "Soergel"
        distance = sum(abs(P-Q))/sum(max(P,Q));
        similarity = 1 - distance;

    case "Kulczynski"
        distance = sum(abs(P-Q))/sum(min(P,Q));
        similarity = 1 / distance;

    case "Canberra"
        distance = sum((abs(P-Q))./(P+Q));
        similarity = nan;

    case "Lorentzian"
        distance = sum(log(1+abs(P-Q)));
        similarity = nan;
    %----------------------------------------------------------------------    

    %----------------------------------------------------------------------
    % Intersection family
    case 'Intersection'
        similarity = sum(min(P,Q));
        distance   = 1 - similarity;

    case 'WaveHedges'
        distance = sum(abs(P-Q)/max(P,Q));
        similarity = 1 - distance;

    case 'Czekanowski'
        similarity = (2*sum(min(P,Q)))/sum(p+Q);
        distance = 1 - similarity;

    case 'Motyka'
        similarity = sum(min(P,Q))/sum(P+Q);
        distance = 1 - similarity;

    case 'Kulczynski'
        similarity = sum(min(P,Q))/sum(abs(P-Q));
        distance = 1 / similarity;

    case 'Ruzicka'
        similarity = sum(min(P,Q))/sum(max(P,Q));
        distance = 1 - similarity;

    case 'Tanimoto'
        distance = sum(max(P,Q)-min(P,Q))/sum(max(P,Q));
        similarity = 1 - distance;
    %----------------------------------------------------------------------   

    %----------------------------------------------------------------------
    % Inner Product family
    case "InnerProduct"
        similarity = sum(P .* Q);
        distance = 1 - similarity;

    case "HarmonicMean"
        similarity = 2 * sum((P .* Q)./(P+Q));
        distance = 1 - similarity;

    case "Cosine"
        similarity = sum(P .* Q) / (sqrt(sum(P.^2))) * (sqrt(sum(Q.^2)));
        distance = 1 - similarity;

    case "Kumar-Hassebrook"
        similarity = sum(P .* Q) / (sum(P.^2) + sum(Q.^2) - sum(P .* Q));
        distance = 1 - similarity;

    case "Jaccard"
        similarity = sum(P .* Q) / (sum(P.^2) + sum(Q.^2) - sum(P .* Q));
        distance = 1 - similarity;

    case "Dice"
        similarity = 2 * sum(P .* Q) / (sum(P.^2) + sum(Q.^2));
        distance   = 1 - similarity;
    %----------------------------------------------------------------------    

    %----------------------------------------------------------------------
    % Fidelity family or Squared-chord family
    case 'Fidelity'
        similarity = sum(sqrt(P .* Q));
        distance = 1 - similarity;

    case 'Bhattacharya'
        distance   = -log(sum(sqrt(P .* Q)));
        similarity = 1 - distance;

    case 'Hellinger'
        distance = 2 * sqrt(1 - sum(sqrt(P .* Q)));
        similarity = nan;

    case 'Matusita'
        distance = sqrt(2 - 2 * sum(sqrt(P .* Q)));
        similarity = 1 - distance;

    case 'Squared-chord'
        distance = sum((sqrt(P)-sqrt(Q)).^2);
        similarity = 1 - distance;
    %----------------------------------------------------------------------    

    %----------------------------------------------------------------------
    % Squared L2 family or χ2 family
    case 'SquaredEuclidean'
        distance = sum((P-Q).^2);
        similarity = 1 - distance;

    case 'Pearson'
        distance = sum(((P-Q).^2./Q));
        similarity = nan;

    case 'Neyman'
        distance = sum(((P-Q).^2./P));
        similarity = nan;

    case 'Squared'
        distance = sum(((P-Q).^2./(P+Q)));
        similarity = 1 - distance;

    case 'ProbabilisticSymmetrix'
        distance = 2*sum(((P-Q).^2./(P+Q)));
        similarity = nan;

    case 'Divergence'
        distance = 2 * sum(((P-Q).^2./(P+Q).^2));
        similarity = nan;

    case 'Clark'
        distance = sqrt(sum((abs(P-Q)./(P+Q)).^2));
        similarity = nan;

    case 'AdditiveSymmetric'
        distance = sum(((P-Q).^2).*(P+Q)./(P .* Q));
        similarity = nan;
    %----------------------------------------------------------------------


    %----------------------------------------------------------------------
    % Shannon’s entropy family    
    case 'Kullback-Leibler'
        distance = sum(P .* log(P./Q));
        similarity = nan;

    case 'Jeffreys'
        distance = sum((P-Q).*log(P./Q));
        similarity = nan;

    case 'KDivergence'
        distance = sum(P .* log((2.*P)./(P+Q)));
        similarity = 1 - distance;

    case 'Topsoe'
        distance = sum(P .* log((2.*P)./(P+Q)) + Q .* log((2.*Q)./(P+Q)));
        similarity = nan;

    case 'Jensen-Shannon'
        distance = (1/2) * (sum(P .* log((2.*P)./(P+Q))) + sum(Q .* log((2.*Q)./(P+Q))));
        similarity = nan;

    case 'JensenDifference'
        distance = sum((P.*log(P) + Q.*log(Q))/2 - ((P+Q)/2 .* log((P+Q)/2)));
        similarity = nan;
    %----------------------------------------------------------------------


    %----------------------------------------------------------------------
    % Combinations    
    case 'Taneja'
        distance = sum(((P+Q)/2) .* (log((P+Q)/(2 * sqrt(P .* Q)))));
        similarity = 1 - distance;

    case 'Kumar-Johnson'
        distance = sum(((P.^2 - Q.^2).^2) ./ (2.*(P .* Q).^(3/2)));
        similarity = nan;

    case 'Average'
        distance = (sum(abs(P-Q)) + max(abs(P-Q))) / 2;
        similarity = 1 - distance;
    %----------------------------------------------------------------------

    otherwise
        error('Invalid Method Name.')
end

% Correlation coefficient of the two distributions
R = corrcoef(P, Q);
corr_coeffient = R(1,2);

end

