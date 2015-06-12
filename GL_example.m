%
%   Multi-Atlas Segmentation as a Graph Labelling Problem: 
%   Application to Partially Annotated Atlas Data
%
%   Example application on a dataset of synthetic tubular structures.
%   This code allows explorations of different graph configurations as
%   proposed in [1].
%   
%
%   Please cite the following paper when using this software:
%   [1] Koch, Lisa M.; Rajchl, M.; Tong, T.; Passerat-Palmbach, J.;
%       Aljabar, P.; Rueckert, D.
%       Multi-Atlas Segmentation as a Graph Labelling Problem: 
%       Application to Partially Annotated Atlas Data 
%       IPMI, 2015
% 
% 
%   Created by lkoch, 2015-06-08
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% your inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Choose graph configuration / experiment type 
%    These are explained in more detail in [1].
type = 'GLa'; % options: 'GLa', 'GLb', 'MAS', 'MASr'

% 2. Choose number of atlases (max. 10)
noAtlases = 5;

% 3. Choose proportion of labelled voxels (recommended >= 0.2)
ratioLabelled = 0.2;

% 4. Specify dataset (example set provided)
fData    = 'synthDataNonlinear.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% graph configurations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. data term: which part of the data is annotated?
%    dataFLAG 1: random slices
%             2: whole images
%             4: evenly spaced slices (with random offset)

% 2. propagation term: how are the images interconnected
%    propagationFLAG 1: all images (target and atlases) are
%                         interconnected
%                    3: all atlases are connected to the target

% 2. regularisation term: which images are regularised
%    regularisationFLAG 0: no regularisation
%                       1: regularisation in target
%                       2: regularisation in all images


if strcmp( type, 'GLa')
    
    % GLa) missing slices, regularisation in all images, propagation
    % between target and atlases
    dataFLAG            = 4;
    propagationFLAG     = 3;
    regularisationFLAG  = 2;

elseif strcmp( type, 'GLb')

    % GLb) missing slices, regularisation in all images, propagation
    % between all images
    dataFLAG            = 4;
    propagationFLAG     = 1;
    regularisationFLAG  = 2;

elseif strcmp( type, 'MAS')

    % MAS) missing images, no regularisation, propagation
    % between target and atlases
    dataFLAG            = 2;
    propagationFLAG     = 3;
    regularisationFLAG  = 0;

elseif strcmp( type, 'MASr')

    % MASr) missing images, regularisation in target, propagation
    % between target and atlases
    dataFLAG            = 2;
    propagationFLAG     = 3;
    regularisationFLAG  = 1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% parameters: adapt if you wish, but not necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% experiment parameters
params = {};
params.visFLAG          = 1;                        % 0=off, 1=on
params.betaType         = 'MeanSquareDiffs';        % propagation weight measure
params.betaScale        = 'local';                  % propagation weight scale
params.mapType          = 'NegativeExpGradient';    % regularisation measure
params.rng_labelselection = 0;

params.numGraphs        = 1 + noAtlases;
params.ratioLabelled    = ratioLabelled;
params.smoothFLAG       = regularisationFLAG;
params.graphconfigFLAG  = propagationFLAG;
params.partialFLAG      = dataFLAG;

params.smoothScale      = 0.3;
params.smoothSigma      = .5;

% convergence parameters
converg = {};
converg.cc              = 0.3;
converg.errbound        = 1e-6;
converg.numIter     	= 50;
converg.steps           = 0.11;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% image and label data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
data = load( fData );
ur   = data.images(:,:,:,1:params.numGraphs);
l    = data.labels(:,:,:,1:params.numGraphs);
gt   = l;

% remove partial labels
l = createPartialLabels( l, params );

figure();
subplot(121); 
imshow(ur(:,:,ceil(end/2),1),[0 max(ur(:))]);
title('target image')

subplot(122); 
imshow(gt(:,:,ceil(end/2),1),[0 1]);
title('ground truth segmentation')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run graph labelling and evaluate results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = graphLabelling( ur, l, params, converg );

figure();
subplot(1,3,1)
p = patch(isosurface(u(:,:,:,1)>0.5));
isonormals(u(:,:,:,1)>0.5,p)
set(p,'FaceColor','red','EdgeColor','none');
title('segmentation result')

daspect([1 1 1])
view(3); axis tight
camlight
lighting gouraud

subplot(1,3,2)
p = patch(isosurface(gt(:,:,:,1)>0));
isonormals(gt(:,:,:,1)>0,p)
set(p,'FaceColor','green','EdgeColor','none');
title('ground truth')

daspect([1 1 1])
view(3); axis tight
camlight
lighting gouraud

subplot(1,3,3)
p = patch(isosurface(u(:,:,:,1)>0.5));
isonormals(u(:,:,:,1)>0.5,p)
set(p,'FaceColor','red','EdgeColor','none');
p = patch(isosurface(gt(:,:,:,1)>0));
isonormals(gt(:,:,:,1)>0,p)
set(p,'FaceColor','green','EdgeColor','none');
title('overlaid')

daspect([1 1 1])
view(3); axis tight
camlight
lighting gouraud

