function [u, erriter, i, timet] = graphLabelling( ur, l, params, converg )
% 
% Performs graph construction (according to parameters) and subsequent 
% optimisation to achieve an optimal labelling on a stack of partially
% labelled images.
% 
% IN:   
%       ur      stack of images
%       l       corresponding stack of labels (partial, with -1 denoting unlabelled
%               voxels)
%       params  experiment parameters
%       converg convergence parameters
% 
% OUT:   
%       u       stack of continuous label maps
%       erriter error
%       i       number of iterations
%       timet   time
% 
%   Please cite the relevant literature when using this software:
%   [1] Koch, Lisa M.; Rajchl, M.; Tong, T.; Passerat-Palmbach, J.;
%       Aljabar, P.; Rueckert, D.
%       Multi-Atlas Segmentation as a Graph Labelling Problem: 
%       Application to Partially Annotated Atlas Data 
%       IPMI, 2015
% 
%   Created by lkoch, 2015-01-28
%   

% Input validation
if nargin < 4
    error('Function requires at least 4 inputs.');
end

%TODO: thorough params/converg validation

algParams = struct( ...,
               'visFLAG', 0, ...
            'smoothFLAG', 1, ...
              'betaType', 'MeanSquareDiffs', ...
             'betaScale', 'local', ...
               'mapType', 'NegativeExpGradient', ...
               'gpuFLAG', 0, ...
             'numGraphs', 11, ...
           'smoothScale', 0.1, ...
           'smoothSigma', 50, ...
       'graphconfigFLAG', 3, ...
               'spacing', [1 1 1] ...
       );

fieldNames = fieldnames(algParams);
paramsNames = fieldnames(params);

for ii = 1:numel(paramsNames)
    
    if any(strcmp(paramsNames{ii},fieldNames))
        % overwrite options. If you want you can test for the right class here
        % Also, if you find out that there is an option you keep getting wrong,
        % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
        algParams.(paramsNames{ii}) = params.(paramsNames{ii});
    else
        warning('%s is not a recognized/necessary parameter name and is ignored', ...
            paramsNames{ii});
    end
end

convParams = struct( ...,
          'cc', 0.3000, ...
    'errbound', 1.0000e-04, ...
     'numIter', 300, ...
       'steps', 0.1100 ...
       );

fieldNames = fieldnames(convParams);
paramsNames = fieldnames(converg);

for ii = 1:numel(paramsNames)
    
    if any(strcmp(paramsNames{ii},fieldNames))
        % overwrite options. If you want you can test for the right class here
        % Also, if you find out that there is an option you keep getting wrong,
        % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
        convParams.(paramsNames{ii}) = converg.(paramsNames{ii});
    else
        error('%s is not a recognized parameter name',paramsNames{ii})
    end
end


%% Run optimiser

[Cs, Ct, beta, alpha]   = graphConstruction(ur, l, algParams);
algParams.uf            = zeros(size(Cs));
algParams.ub            = ones(size(Cs));
algParams.Cs            = Cs;
algParams.Ct            = Ct;
algParams.beta          = beta;
algParams.alpha         = alpha;

[u, erriter, i, timet]  = graphOptimisation(algParams, convParams);

end