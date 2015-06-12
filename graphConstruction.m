function [ Cs, Ct, beta, alpha ] = graphConstruction( ur, l, params )
%
% Performs graph construction (data term, propagation term, regularisation
% term)
%	
% IN:   
%       ur      stack of images
%       l       corresponding stack of labels (partial, with -1 denoting unlabelled
%               voxels)
%       params  experiment parameters
% 
% OUT:   
%       Cs, Ct  source and sink connections (data term)
%       beta    inter-image connections (propagation term)
%       alpha   regularisation constraints
% 
%   Created by lkoch, 2015-01-28
%   

%% input validation
if nargin < 3
    error('Function requires at least three inputs.');
end

% regularisation parameters
if ~isfield(params, 'smoothScale')
    params.smoothScale = .1;
end
if ~isfield(params, 'smoothSigma')
    params.smoothSigma = 200;
end

%% sink and source flows for atlases

% set labels 1 and 0 (fg and bg) as terminal connections
% unlabelled voxels should have been labelled as -1 and will then remain
% unconnected
Cs = (l==1) * 100;
Ct = (l==0) * 100;

%% inter-image connections

[rows, cols, heights] = size(Cs(:,:,:,1));
beta            = zeros(rows,cols,heights,params.numGraphs,params.numGraphs);

if params.graphconfigFLAG == 1
    
    for ii=1:params.numGraphs
        for jj=1:params.numGraphs
            
            if ii==jj
                continue
            end
            
            beta(:,:,:,ii,jj)   = fusionWeights( ur(:,:,:,ii), ur(:,:,:,jj), ...
                'SimilarityMeasure', params.betaType, 'Scale', params.betaScale );
            
        end
    end
    
elseif params.graphconfigFLAG == 3
    
    for jj=2:params.numGraphs                
        beta(:,:,:,1,jj) = fusionWeights( ur(:,:,:,1), ur(:,:,:,jj), ...
            'SimilarityMeasure', params.betaType, 'Scale', params.betaScale );
        
        beta(:,:,:,1,jj) = fusionWeights( ur(:,:,:,1), ur(:,:,:,jj), ...
            'SimilarityMeasure', params.betaType, 'Scale', params.betaScale );
        beta(:,:,:,jj,1) = beta(:,:,:,1,jj);
    end
    
end

%% intra-image connections

alpha = zeros(rows,cols,heights,params.numGraphs);

if params.smoothFLAG==1
    % smoothing in target only
    alpha(:,:,:,1) = smoothingWeights(ur(:,:,:,1), ...
        'MapType', params.mapType, ...
        'ScalingCoefficient', params.smoothScale, ...
        'Sigma', params.smoothSigma);
    
elseif params.smoothFLAG==2
    % smoothing in all images
    for jj=1:params.numGraphs
        alpha(:,:,:,jj) = smoothingWeights(ur(:,:,:,jj), ...
            'MapType', params.mapType, ...
            'ScalingCoefficient', params.smoothScale, ...
            'Sigma', params.smoothSigma);        
    end   
end

end


