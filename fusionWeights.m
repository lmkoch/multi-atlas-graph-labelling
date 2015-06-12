function [ beta ] = fusionWeights( img1, img2, varargin )
%fusionWeights calculates similarity weights between two images
%   
% IN:   
%       img1    images to compare
%       img2    
% 
% OUT:   
%       beta    weights between each pair of voxels in img1 and img2
%               same dimension as img1 and img2
%	
%   Created by lkoch, 2014-11-08
%   

% Input validation
if nargin < 2
    error('Function requires at least two inputs.');
end

assert( isequal(size(img1),size(img2)), 'img1 and img2 should be of equal dimensions.' );

options = struct( ...
    'RegionOfInterest', ones(size(img1))>0, ...
    'Scale', 'global', ...
    'SimilarityMeasure', 'CrossCorrelation' ...
    );

optionNames = fieldnames(options);

nArgs = length(varargin);
assert( round(nArgs/2)==nArgs/2, 'weightMaps needs propertyName/propertyValue pairs');

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    
    if any(strcmp(pair{1},'SimilarityMeasure'))
        options.(pair{1}) = validatestring(pair{2}, ...
            {'CrossCorrelation','SumOfSquareDiffs','MeanSquareDiffs','Bai2013','Uniform'} ...
            );
        
    elseif any(strcmp(pair{1},'Scale'))
        options.(pair{1}) = validatestring(pair{2},{'global','local'});
        
    elseif any(strcmp(pair{1},'RegionOfInterest'))
        options.(pair{1}) = pair{2}>0;
        
    elseif any(strcmp(pair{1},optionNames))
        % overwrite options. If you want you can test for the right class here
        % Also, if you find out that there is an option you keep getting wrong,
        % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name',pair{1})
    end
end


% Calculate weight  maps

beta = zeros(size(img1));

switch options.Scale
    case 'global'
        
        switch options.SimilarityMeasure
            case 'CrossCorrelation'
                roi = options.RegionOfInterest;
                beta(roi) = corr(img1(roi), img2(roi)) * (corr(img1(roi), img2(roi))>0);
                
            case 'SumOfSquareDiffs'
                roi = options.RegionOfInterest;
                beta(roi) = (img1(roi)-img2(roi)).^2;
                
            case 'Uniform'
                beta(:) = 1;
        end
        
    case 'local'
        
        switch options.SimilarityMeasure
            case 'SumOfSquareDiffs'
                roi = options.RegionOfInterest;
                
                patchsize = 9;
                patch = ones( patchsize * ones(1,ndims(img1)) );
                patch = patch / sum(patch(:));
                
                tmp = convn( (img1-img2).^2, patch, 'same');
                beta(roi) = tmp(roi).^-.2;
                
            case 'MeanSquareDiffs'
                % Parameters according to Artaechevarria et al., TMI 2009
                % Combination strategies in multi-atlas image segmentation: 
                % Application to brain MR data
                
                roi = options.RegionOfInterest;
                
                patchsize   = 9;
%                 gain        = -1;
                gain        = -.2;
                patch       = ones( patchsize * ones(1,ndims(img1)) );
                patch       = patch / sum(patch(:));
                
                tmp         = convn( (img1-img2).^2, patch, 'same');
                beta(roi)   = tmp(roi).^gain;
                
            case 'Bai2013'
                % Parameters according to Bai et al., TMI 2013
                % A Probabilistic Patch-Based Label Fusion Model
                % for Multi-Atlas Segmentation With Registration
                % Refinement: Application to Cardiac MR Images
                
                roi = options.RegionOfInterest;
                
                sigma1      = 50;                
                patch       = ones( [3 3 1] );                
                
                tmp         = convn( (img1-img2).^2, patch, 'same');
                beta(roi)   = 1/(sqrt(2*pi) * sigma1) * ...
                    exp( -tmp(roi) / ( 2*sigma1^2*sum(patch(:))) );
                
        end
        
end

