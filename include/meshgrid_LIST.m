function [xx,yy,zz,xxyyzz] = meshgrid_LIST(x,y,z)
zz=[];xx=[];yy=[];
wxxyyzz=[];
if nargin == 0 % || nargin > 1%(nargin > 1 && nargout > nargin)
    error(message('MATLAB:meshgrid:NotEnoughInputs'));
end

if nargin == 2 || (nargin == 1 && nargout < 3) % 2-D array case
    if nargin == 1
        y = x;
    end
    if isempty(x) || isempty(y)
        xx = zeros(0,class(x));
        yy = zeros(0,class(y));
    else
        xrow = full(x(:)).'; % Make sure x is a full row vector.
        ycol = full(y(:));   % Make sure y is a full column vector.
        xx = repmat(xrow,size(ycol));
        yy = repmat(ycol,size(xrow));
    end
else  % 3-D array case
    if nargin == 1
        y = x;
        z = x;
    end
    if isempty(x) || isempty(y) || isempty(z)
        xx = zeros(0,class(x));
        yy = zeros(0,class(y));
        zz = zeros(0,class(z));
    else
        nx = numel(x);
        ny = numel(y);
        nz = numel(z);
        xx = reshape(full(x),[1 nx 1]); % Make sure x is a full row vector.
        yy = reshape(full(y),[ny 1 1]); % Make sure y is a full column vector.
        zz = reshape(full(z),[1 1 nz]); % Make sure z is a full page vector.
        xx = repmat(xx, ny, 1, nz);
        yy = repmat(yy, 1, nx, nz);
        zz = repmat(zz, ny, nx, 1);
    end
end
try
xx=reshape(xx,[],1);end
try
yy=reshape(yy,[],1);
end
try
zz=reshape(zz,[],1);
end
try 
xxyyzz=xx;
end
try
    xxyyzz=[xx,yy];end
try
    xxyyzz=[xx,yy,zz];end

