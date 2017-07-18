function [Tx, Ty]=demosimp(im, refim)
% Original image size 
dimen=size(refim);

% Go across sublevels
Tx=0;Ty=0;
for level=3:-1:1
    
    % Size at the smallest hierarchical level, when we resize to smaller
    M=floor(dimen/2^(level-1));

    % current image size
    refim1 = imresize(refim,M); 
    im1 = imresize(im,M);
    im2 = movepixels(im1,Tx,Ty); 
    
    [Tx1,Ty1] = demos(im2,refim1);  
    Tx=Tx+Tx1;
    Ty=Ty+Ty1;
    if level>1
        M=floor(dimen/2^(level-2));
        Tx=imresize(Tx,M);
        Ty=imresize(Ty,M);
    end
end