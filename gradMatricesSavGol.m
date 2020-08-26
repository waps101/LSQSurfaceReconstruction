function [Dx,Dy,Serr,S] = gradMatricesSavGol(mask,nOrder,nSize,verbose)
%GRADMATRICESSAVGOL Build numerical derivative and smoothness matrices
%   For an arbitrary foreground mask, build derivative and smoothness prior
%   matrices based on 2D Savitzky-Golay filters of order nOrder and size
%   nSize x nSize. Pixels on or near the boundary where the filter would
%   cover non-foreground pixels use an irregularly shaped filter comprising
%   the nSize^2 closest foreground pixels.
%
% Inputs:
%    mask   - H x W logical mask indicating which pixels are foreground
%    nOrder - order of 2D polynomial to use
%    nSize  - square filter dimension for non-boundary pixels
% Outputs:
%    Dx     - sparse matrix of size npix x npix where npix is the number of
%             foreground pixels evaluates the x derivative when multiplied
%             by the vectorised version of the function values
%    Dy     - same as Dx for y derivative
%    Serr   - same size as Dx but measures deviation from the smoothed 
%             version of the function, i.e. a version which is locally LSQ 
%             fitted by a polynomial
%    S      - same size as Dx but computes smoothed version of the 
%             function, i.e. a version which is locally LSQ fitted by a 
%             polynomial
%
% William Smith
% University of York
% 2018

if nargin<4
    verbose = false;
end

if verbose
    disp(['Using ' num2str(nSize) 'x' num2str(nSize) ' kernel with order ' num2str(nOrder) ' polynomial surface']);
end

rows = size(mask,1);
cols = size(mask,2);

% The number of usable pixels
npix = sum(mask(:));

% Build lookup table relating x,y coordinate of valid pixels to index
% position in vectorised representation
indices = zeros(size(mask));
indices(mask) = 1:npix;

% Precompute the standard Savitzky Golay filters used for internal pixels
Terms = [];
XPow = [];
YPow = [];
syms x y real;
for j=0:nOrder
    for i=0:nOrder-j
        Terms = [Terms; (x^i)*(y^j)];
        XPow = [XPow; i];
        YPow = [YPow; j];
    end
end
% Compute A matrix for a nSize x nSize window
A = [];
for x = -(nSize-1)/2:(nSize-1)/2 % important to loop through in same scan order as image patch pixels
    for y = (nSize-1)/2:-1:-(nSize-1)/2
        A = [A; subs(Terms')];
    end
end
% Compute coefficient matrix
C = inv(A'*A)*A';
% Pull out coefficients
SG = [];
nTerms = size(Terms,1);
for i=1:nTerms
    SG(:,:,i) = reshape(C(i,:),[nSize,nSize]);
end

dx = squeeze(SG(:,:,XPow==1 & YPow==0));
dy = squeeze(SG(:,:,XPow==0 & YPow==1));
smooth = squeeze(SG(:,:,XPow==0 & YPow==0));
% Subtract 1 from the weight of pixel itself so that the matrix measures 
% deviation from the smoothed value rather than giving the smoothed value
% itself:
smooth((nSize+1)/2,(nSize+1)/2)=smooth((nSize+1)/2,(nSize+1)/2)-1;

% Find non-boundary pixels where entire mask can be centered over pixel
internal = imfilter(uint16(mask),ones(nSize,nSize))==nSize^2;

% Initialise and preallocate sparse for sparse entries
NumEq=0;

ix = zeros(npix*nSize^2,1);
jx = zeros(npix*nSize^2,1);
sx = zeros(npix*nSize^2,1);
kx=0;

iy = zeros(npix*nSize^2,1);
jy = zeros(npix*nSize^2,1);
sy = zeros(npix*nSize^2,1);
ky=0;

is = zeros(npix*nSize^2,1);
js = zeros(npix*nSize^2,1);
ss = zeros(npix*nSize^2,1);
ks=0;

[maskr,maskc]=find(mask);
% Put mask pixels into a KD tree so that NN searches are fast later
ns = createns([maskr maskc],'nsmethod','kdtree');

if verbose
    disp([num2str(npix) ' foreground pixels']);
    disp([num2str(sum(internal(:))) ' have full neighbourhood, ' num2str(npix-sum(internal(:))) ' are boundary cases']);
end

% This is much faster than using subs inside the loop!
Terms = matlabFunction(Terms');

for col=1:cols
    for row=1:rows
        if internal(row,col)
            % Can use standard filters
            NumEq = NumEq+1;
            for c=-(nSize-1)/2:(nSize-1)/2
                for r=-(nSize-1)/2:(nSize-1)/2                    
                    % Add equation to the x derivative matrix
                    if dx(r+(nSize-1)/2+1,c+(nSize-1)/2+1)~=0
                        kx = kx+1;
                        ix(kx)=NumEq; jx(kx)=indices(row+r,col+c); sx(kx)=dx(r+(nSize-1)/2+1,c+(nSize-1)/2+1);
                    end
                    % Add equation to the y derivative matrix
                    if dy(r+(nSize-1)/2+1,c+(nSize-1)/2+1)~=0
                        ky = ky+1;
                        iy(ky)=NumEq; jy(ky)=indices(row+r,col+c); sy(ky)=dy(r+(nSize-1)/2+1,c+(nSize-1)/2+1);
                    end
                    % Add equation to the smoothness matrix
                    if smooth(r+(nSize-1)/2+1,c+(nSize-1)/2+1)~=0
                        ks = ks+1;
                        is(ks)=NumEq; js(ks)=indices(row+r,col+c); ss(ks)=smooth(r+(nSize-1)/2+1,c+(nSize-1)/2+1);
                    end
                end
            end
        elseif mask(row,col)
            % Pixel is inside the mask but we can't use the standard
            % filters - build a custom one
            
            % 1. Find the nSize^2 closest pixels (including self)
            idx = knnsearch(ns,[row col],'K',nSize^2);
            % break ties by taking first nSize^2 neighbours
            idx = idx(1:nSize^2);
            % 2. Build the custom Savitzky-Golay filters for this pixel in
            % a coordinate system where this pixel is (0,0)
            A = [];
            for i=1:nSize^2
                x = maskc(idx(i))-col;
                y = row-maskr(idx(i));
                A = [A; Terms(x,y)];
            end
            % Compute coefficient matrix
            C = inv(A'*A)*A';
            NumEq = NumEq+1;
            % Subtract 1 from the weight of pixel itself so that the matrix
            % measures deviation from the smoothed value rather than giving
            % the smoothed value itself
            C(XPow==0 & YPow==0,1) = C(XPow==0 & YPow==0,1)-1;
            % 3. Add the equations for this filter to the sparse matrices
            for i=1:nSize^2
                % Add equation to the x derivative matrix
                if C(XPow==1 & YPow==0,i)~=0
                    kx = kx+1;
                    ix(kx)=NumEq; jx(kx)=indices(maskr(idx(i)),maskc(idx(i))); sx(kx)=C(XPow==1 & YPow==0,i);
                end
                % Add equation to the y derivative matrix
                if C(XPow==0 & YPow==1,i)~=0
                    ky = ky+1;
                    iy(ky)=NumEq; jy(ky)=indices(maskr(idx(i)),maskc(idx(i))); sy(ky)=C(XPow==0 & YPow==1,i);
                end
                % Add equation to the smoothness matrix
                if C(XPow==0 & YPow==0,i)~=0
                    ks = ks+1;
                    is(ks)=NumEq; js(ks)=indices(maskr(idx(i)),maskc(idx(i))); ss(ks)=C(XPow==0 & YPow==0,i);
                end
            end

        end
    end
end

% Throw away unused elements of the sparse indices
ix=ix(1:kx,1); jx=jx(1:kx,1); sx=sx(1:kx,1);
iy=iy(1:ky,1); jy=jy(1:ky,1); sy=sy(1:ky,1);
is=is(1:ks,1); js=js(1:ks,1); ss=ss(1:ks,1);

% Build the sparse matrices
Dx = sparse(ix,jx,sx,npix,npix);
Dy = sparse(iy,jy,sy,npix,npix);
S  = sparse(is,js,ss,npix,npix);
Serr = S;
S = S+speye(npix);

end