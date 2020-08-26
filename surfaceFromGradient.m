function [z] = surfaceFromGradient(Dx,Dy,p,q,mask,options)
%HEIGHTFROMGRADIENT Perform least squares surface integration
%   Detailed explanation goes here

% Setup some default options
if (nargin<6) || ~isfield(options,'cameramodel')
    options.cameramodel = 'orthographic';
end
if ~isfield(options,'precond')
    options.precond = 'none';
end
if ~isfield(options,'solver')
    options.solver = 'direct';
end
if ~isfield(options,'structure')
    options.structure = 'stacked';
end
if ~isfield(options,'regularisation')
    options.regularisation = 'none';
end
if ~isfield(options,'z0')
    options.z0 = zeros(size(mask));
end

if strcmp(options.structure,'stacked')
    if strcmp(options.solver,'pcg') || strcmp(options.precond,'ichol')
        error('Can only use PCG solver or ichol preconditioning with a square matrix');
    end
end

if strcmp(options.cameramodel,'orthographic')
    npix = sum(mask(:));
    norms = sqrt(p.^2+q.^2+1);
    Nx = -p./norms;
    Ny = -q./norms;
    Nz = 1./norms;
    if strcmp(options.structure,'stacked')
        A = [Dx;Dy];
        %A = [sparse(1:npix,1:npix,Nz(mask),npix,npix)*Dx;sparse(1:npix,1:npix,Nz(mask),npix,npix)*Dy];
        b = [p(mask);q(mask)];
        %b = [-Nx(mask);-Ny(mask)];
    elseif strcmp(options.structure,'square')
        A = Dx'*Dx + Dy'*Dy;
        b = Dx'*p(mask) + Dy'*q(mask);
    end
elseif strcmp(options.cameramodel,'perspective')
    [x,y]=meshgrid(1:size(mask,2),size(mask,1):-1:1);
    npix = sum(mask(:));
    X = sparse(1:npix,1:npix,x(mask)-options.cx,npix,npix);
    Y = sparse(1:npix,1:npix,y(mask)-options.cy,npix,npix);
    % Build tangent vector matrices
    Tx = [-1/options.fx.*X -1/options.fx.*speye(npix); ...
        -1/options.fy.*Y sparse([],[],[],npix,npix); ...
        speye(npix) sparse([],[],[],npix,npix)] * [Dx; speye(npix)];
    Ty = [-1/options.fx.*X sparse([],[],[],npix,npix); ...
        -1/options.fy.*Y -1/options.fy.*speye(npix); ...
        speye(npix) sparse([],[],[],npix,npix)] * [Dy; speye(npix)];
    % Build matrix for performing dot product with target surface normals
    norms = sqrt(p.^2+q.^2+1);
    Nx = -p./norms;
    Ny = -q./norms;
    Nz = 1./norms;
    N = sparse([1:npix 1:npix 1:npix],1:3*npix,[Nx(mask); Ny(mask); Nz(mask)],npix,3*npix);
    
    if strcmp(options.structure,'stacked')
        A = [N*Tx; N*Ty];
        b = zeros(2*npix,1);
    elseif strcmp(options.structure,'square')
        A = (N*Tx)'*(N*Tx) + (N*Ty)'*N*Ty;
        b = zeros(npix,1);
    end
end

if strcmp(options.regularisation,'z0')
    % options.z0 contains a guide depth map
    if strcmp(options.structure,'stacked')
        A = [A;options.lambda.*speye(size(Dx,1))];
        b = [b;options.z0(mask)];
    elseif strcmp(options.structure,'square')
        A = A + options.lambda^2.*speye(size(Dx,1));
    end
elseif strcmp(options.regularisation,'smoothness')
    % options.S contains a matrix that measures deviation from smooth
    % surface and options.lambda is the regularisation weight
    if strcmp(options.structure,'stacked')
        A = [A;options.lambda.*options.S];
        b = [b;zeros(size(options.S,1),1)];
    elseif strcmp(options.structure,'square')
        A = A + options.lambda^2.*(options.S'*options.S);
    end
end

if ~strcmp(options.regularisation,'z0')
    % Need to resolve constant of integration ambiguity
    if strcmp(options.structure,'stacked')
        A = [A;sparse(1,1,1,1,size(Dx,2))];
        b = [b;1];
    elseif strcmp(options.structure,'square')
        A = A + sparse(1,1,1,size(Dx,1),size(Dx,2));
    end
end

if(strcmp(options.precond,'none'))
    precondL = [];
    precondR = [];
elseif strcmp(options.precond,'ichol')
    precondL = ichol(A,struct('type','ict','droptol',1e-03,'michol','on'));
    precondR = precondL';
end
% Resolution
z = options.z0;
if(strcmp(options.solver,'direct')) % Calls cholesky
    z(mask) = A\b;
elseif(strcmp(options.solver,'pcg')) % Calls CG
    z(mask) = pcg(A,b,1e-4,1000,precondL,precondR,z(mask));
end
z(~mask) = NaN;

end

