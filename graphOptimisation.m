function [u, erriter, i, timet] = graphOptimisation(params, converg)
%graphOptimisation continuous max-flow (CMF) based optimisation
%
% IN:   
%       params  experiment parameters
%       converg convergence parameters
% 
% OUT:   
%       u       stack of continuous label maps
%       erriter error
%       i       number of iterations
%       timet   time
% 
% 
%   Please cite the relevant literature when using this software.
%   This version of the CMF was proposed in the following paper:
%   [1] Koch, Lisa M.; Rajchl, M.; Tong, T.; Passerat-Palmbach, J.;
%       Aljabar, P.; Rueckert, D.
%       Multi-Atlas Segmentation as a Graph Labelling Problem: 
%       Application to Partially Annotated Atlas Data 
%       IPMI, 2015
% 
%   The original algorithm was proposed in the following paper:
%   [2] Yuan, J.; Bae, E.;  Tai, X.-C. 
%       A Study on Continuous Max-Flow and Min-Cut Approaches 
%       CVPR, 2010
% 
%   This implementation relies heavily on the implementation of [2]
%   provided here: https://sites.google.com/site/wwwjingyuan/
% 
% 
%   Created by lkoch, 2015-01-28
%


% experiment parameters
Cs          = params.Cs;
Ct          = params.Ct;
gpuFLAG     = params.gpuFLAG;
numGraphs   = params.numGraphs;
beta        = params.beta;
alpha       = params.alpha;
visFLAG     = params.visFLAG;


% convergence parameters
cc          = converg.cc;
errbound    = converg.errbound;
numIter     = converg.numIter;
steps       = converg.steps;


[rows, cols, heights] = size(Cs(:,:,:,1));
imgSize = rows*cols*heights;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = ones(size(Cs)) *.5; % this ensures .5 for ambivalent outputs
ps = min(Cs, Ct);
pt = ps;

if gpuFLAG
    
    alpha   = gpuArray(alpha);
    Cs      = gpuArray(Cs);
    Ct      = gpuArray(Ct);
    
    u = gpuArray(u);
    ps = gpuArray(ps);
    pt = gpuArray(pt);
    
    px = gpuArray(zeros(rows, cols+1, heights,numGraphs));
    py = gpuArray(zeros(rows+1, cols, heights,numGraphs));
    pz = gpuArray(zeros(rows, cols, heights+1,numGraphs));
    
    divp = gpuArray(zeros(rows, cols, heights,numGraphs));
    
    erru = gpuArray(zeros(rows,cols,heights,numGraphs));
    erriter = gpuArray(zeros(numIter,1));
    
    beta    = gpuArray(beta);
    
else
    
    px = zeros(rows, cols+1, heights,numGraphs, 'single');
    py = zeros(rows+1, cols, heights,numGraphs, 'single');
    pz = zeros(rows, cols, heights+1,numGraphs, 'single');
    
    divp = zeros(rows, cols, heights,numGraphs, 'single');
    
    erru = zeros(rows,cols,heights,numGraphs);
    erriter = zeros(numIter,1);
    
end

for j=1:numGraphs
    divp(:,:,:,j) = - py(1:rows,:,:,j) + py(2:rows+1,:,:,j) + px(:,2:cols+1,:,j) ...
        - px(:,1:cols,:,j) + pz(:,:,2:heights+1,j) - pz(:,:,1:heights,j);
end

if gpuFLAG
    qq = gpuArray(zeros(rows,cols,heights,numGraphs,numGraphs));
else
    qq = zeros(rows,cols,heights,numGraphs,numGraphs, 'single');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visFLAG
    figure('Units', 'pixels', ...
        'Position', [100 100 900 375]);
end

tOverall=tic;
tIter=tic;
for i = 1:numIter
    
    % For each subgraph (image), update all spatial flows (1.), then update
    % inter-image flows to all neighbouring subgraphs (2.), then update
    % source/sink flows
    
    for j = 1:numGraphs
        
        beta_j = beta(:,:,:,:,j);
        qq_j    = qq(:,:,:,j,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. update spatial flows
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % the spatial flows need to be updated only if the alphas are non-zero..
        % otherwise nothing will happen anyway
        if nnz(alpha(:,:,:,j))>0
            
            pts = divp(:,:,:,j) - ps(:,:,:,j) + pt(:,:,:,j) - u(:,:,:,j)/cc;
            
            % note the sign
            pts = pts + sum(qq_j(:,:,:,1,:), 5);
            
            px(:,2:cols,:,j) = px(:,2:cols,:,j) + steps*(pts(:,2:cols,:) - pts(:,1:cols-1,:));
            py(2:rows,:,:,j) = py(2:rows,:,:,j) + steps*(pts(2:rows,:,:) - pts(1:rows-1,:,:));
            pz(:,:,2:heights,j) = pz(:,:,2:heights,j) + steps*(pts(:,:,2:heights) - pts(:,:,1:heights-1));
            
            % the following steps give the projection to make |p(x)| <= alpha(x)
            
            gk = sqrt((px(:,1:cols,:,j).^2 + px(:,2:cols+1,:,j).^2 + ...
                py(1:rows,:,:,j).^2 + py(2:rows+1,:,:,j).^2 + ...
                pz(:,:,1:heights,j).^2 + pz(:,:,2:heights+1,j).^2)*0.5);
            
            gk = double(gk <= alpha(:,:,:,j)) + (double(~(gk <= alpha(:,:,:,j)))+eps).*((gk+eps) ./ alpha(:,:,:,j));
            gk = 1 ./ gk;
            
            px(:,2:cols,:,j) = (0.5*(gk(:,2:cols,:) + gk(:,1:cols-1,:))).*px(:,2:cols,:,j);
            py(2:rows,:,:,j) = (0.5*(gk(2:rows,:,:) + gk(1:rows-1,:,:))).*py(2:rows,:,:,j);
            pz(:,:,2:heights,j) = (0.5*(gk(:,:,2:heights) + gk(:,:,1:heights-1))).*pz(:,:,2:heights,j);
            
            divp(:,:,:,j) = - py(1:rows,:,:,j) + py(2:rows+1,:,:,j) + px(:,2:cols+1,:,j) ...
                - px(:,1:cols,:,j) + pz(:,:,2:heights+1,j) - pz(:,:,1:heights,j);
            
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. update inter-image flows
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for ii = 1:numGraphs
            
            if ii==j
                continue
            end
            
            % this only needs to be done if the betas are non-zero..
            % otherwise nothing will happen anyway
            if nnz(beta_j(:,:,:,ii))==0
                continue
            end
                        
            qq_i    = qq(:,:,:,ii,:);
            qq_j    = qq(:,:,:,j,:);
            
            % update the relaxed flow qq
            
            con1 = divp(:,:,:,ii) - ps(:,:,:,ii) + pt(:,:,:,ii)  - u(:,:,:,ii)/cc + ...
                sum(qq_i(:,:,:,1,:), 5) - qq_i(:,:,:,1,j) ;
            
            con2 = divp(:,:,:,j) - ps(:,:,:,j) + pt(:,:,:,j)  - u(:,:,:,j)/cc + ...
                sum(qq_j(:,:,:,1,:), 5) - qq_j(:,:,:,1,ii) ;
            
            qq_i(:,:,:,1,j) = (con2-con1)/2;
            
            % bidirectional flow constraint
            qq_i(:,:,:,1,j) = max( -beta_j(:,:,:,ii), min( beta_j(:,:,:,ii), qq_i(:,:,:,1,j)) );
            qq_j(:,:,:,1,ii) = -qq_i(:,:,:,1,j);
                        
            qq(:,:,:,ii,:) = gather(qq_i);
            qq(:,:,:,j,:) = gather(qq_j);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3. update terminal flows
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % update the sink flow pt
        pts = - divp(:,:,:,j) + ps(:,:,:,j) + u(:,:,:,j)/cc - sum(qq_j(:,:,:,1,:), 5);
        pt(:,:,:,j) = min(pts, Ct(:,:,:,j));
        
        % update the source flow ps
        pts = divp(:,:,:,j) + pt(:,:,:,j) - u(:,:,:,j)/cc + sum(qq_j(:,:,:,1,:), 5) + 1/cc;
        ps(:,:,:,j) = min(pts, Cs(:,:,:,j));
        
    end
    
    % update the multiplier u for all subgraphs
    for j= 1:numGraphs
        erru(:,:,:,j) = cc*(divp(:,:,:,j) + pt(:,:,:,j)  - ps(:,:,:,j) + sum(qq(:,:,:,j,:), 5) );
    end
    
    u = u - erru;
    
    % evaluate the average error    
    erriter(i) = sum(sum(sum(sum(abs(erru)))))/imgSize/numGraphs;
    
    fprintf('Ending %d after %.2f minutes, err:%g\n', i, toc(tIter)/60, erriter(i))
    
    if visFLAG        
        if numGraphs>=5
            for k=1:numGraphs
                subplot(1,2,1); 
                imshow(u(:,:,ceil(heights/2),1),[0 1]);
                title('estimated labelling')
            end
            subplot(1,1+1,1+1); 
            loglog(erriter); 
            title('error')
            xlabel('iterations')
            ylabel('error')
            drawnow();
            
        else
            for k=1:numGraphs
                subplot(1,numGraphs+1,k); imshow(u(:,:,ceil(heights/2),k),[0 1]);
                title('estimated labelling')
            end
            subplot(1,numGraphs+1,numGraphs+1); loglog(erriter); drawnow();
            title('error')
        end
    end
    
    if (erriter(i) < errbound)
        break;
    end
end

u=gather(u);
erriter=gather(erriter);

timet = toc(tOverall);


if gpuFLAG
    reset(gpuDevice())
end

end