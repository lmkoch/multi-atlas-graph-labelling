function [ alpha ] = smoothingWeights( img, varargin )
%
% Calculates regularisation constraints per voxel
%   
% IN:   
%       img     image to regularise
% 
% OUT:   
%       alpha   regularisation weight for each voxel in img
%               same dimension as img
%	
%   Created by lkoch, 2014-11-08
%   

% Input validation
if nargin < 1
    error('Function requires at least one input.');
end

options = struct( ...
    'MapType', 'Uniform', ...
    'ScalingCoefficient', 1, ...
    'Sigma', 200, ...
    'VoxelSpacing', [1 1 1] ...
    );

optionNames = fieldnames(options);

nArgs = length(varargin);
assert( round(nArgs/2)==nArgs/2, 'weightMaps needs propertyName/propertyValue pairs');

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    
    if any(strcmp(pair{1},'MapType'))
        options.(pair{1}) = validatestring(pair{2},{'InverseGradient','Uniform', 'NegativeExpGradient'});
                
    elseif any(strcmp(pair{1},'ScalingCoefficient'))
        options.(pair{1}) = pair{2};
                
    elseif any(strcmp(pair{1},'Sigma'))
        options.(pair{1}) = pair{2};
                
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

alpha = zeros(size(img));

switch options.MapType
    case 'Uniform'
        alpha(:) = 1;
        
    case 'InverseGradient'
        [dx,dy,dz] = gradient(img);
        dx = dx / options.VoxelSpacing(1);
        dy = dy / options.VoxelSpacing(2);
        dz = dz / options.VoxelSpacing(3);
        
        Gmag = sqrt( dx.^2 + dy.^2 + dz.^2 );                
        alpha = 1-Gmag/max(Gmag(:)); 
        
    case 'NegativeExpGradient'
        [dx,dy,dz] = gradient(img);
        dx = dx / options.VoxelSpacing(1);
        dy = dy / options.VoxelSpacing(2);
        dz = dz / options.VoxelSpacing(3);
        
        Gmag = sqrt( dx.^2 + dy.^2 + dz.^2 );   
        
        sigma = options.Sigma;
        coeff = options.ScalingCoefficient;
        alpha = coeff * exp(- Gmag.^2/ (2*sigma^2) );
end


end

