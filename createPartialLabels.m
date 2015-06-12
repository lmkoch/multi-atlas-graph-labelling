function [ l ] = createPartialLabels( l, params )
% 
% Simulates partially annotated label maps by removing labels from 
% certain voxels according to a chosen configuration (e.g. missing slices,
% missing images).
% 
% IN:   
%       l       stack of labels (partial, with -1 denoting unlabelled
%               voxels)
%       params  experiment parameters
% 
% OUT:   
%       l       stack of partially annotated labels (-1 denoting unlabelled
%               voxels)
% 
%	
%   Created by lkoch, 2015-01-28
%   


% fixing random number generator
rng(params.rng_labelselection); 

[rows, cols, heights] = size(l(:,:,:,1));

for i=1:params.numGraphs
            
    if params.partialFLAG==1    
        % remove random slices
        
        idx = randperm(rows)>(rows*params.ratioLabelled);
        l(idx,:,:,i) = -1;

    elseif params.partialFLAG==3    
        % remove evenly spaced slices with a regular offset
        
        offset = mod( floor( (i-1)/(params.numGraphs * params.ratioLabelled) ), ...
            round(1/params.ratioLabelled) );
        
        spacing = 1/(1-params.ratioLabelled);
        idx = round(spacing:spacing:rows);
        idx = idx + offset;
        
        offset = offset+1;
        if offset>round(1/params.ratioLabelled)
            offset=0;
        end        
        
        idx = idx(idx<=rows);        
        l(idx,:,:,i) = -1;

    elseif params.partialFLAG==4    
        % remove evenly spaced slices with a random offset
                    
        spacing = 1/(1-params.ratioLabelled);
        idx = round(-spacing:spacing:rows);
        offset = randsample(round(1/params.ratioLabelled),1)-1;
        idx = idx + offset;
        
        idx = idx(idx>0);
        idx = idx(idx<=rows);
        
        l(idx,:,:,i) = -1;

    end   
    
end

% remove labels for target image
if params.graphconfigFLAG == 3 || params.graphconfigFLAG == 1   
    l(:,:,:,1) = -1;
end

if params.partialFLAG==2
    % remove random atlases    
    idx = randperm(params.numGraphs-1) > ((params.numGraphs-1)*params.ratioLabelled);    
    l(:,:,:,find(idx)+1) = -1;
end

if params.visFLAG
    figure();
    for k=1:params.numGraphs
        subplot(1,params.numGraphs,k); 
        imshow(l(:,:,ceil(heights/2),k),[min(unique(l)) max(unique(l))]);
        title( sprintf('label %d', k) );
    end
end

end

