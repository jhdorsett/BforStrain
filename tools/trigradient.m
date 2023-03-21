function [dfx,dfy]=trigradient(tri,x,y,z,opt)
%TRIGRADIENT Triangular Gradient.
% [dFx,dFy] = TRIGRADIENT(TRI,X,Y,Z) or TRIGRADIENT(TRI,X,Y,Z,'vertex')
% returns numerical gradients of the function z = F(x,y), dFx = dF/dx
% and DFy = dF/dy at the vertices of the triangles specified in TRI.
% dFx and dFy are column vectors having the same number of elements as X,
% Y, and Z.
% TRI, X, Y, and Z define a triangulation where the triangles are defined
% by the M-by-3 face matrix TRI, such as that returned by DELAUNAY. Each
% row of TRI contains indices into the X, Y, and Z vertex vectors to define
% a single triangular face.
% The gradients are computed by an inverse distance method, whereby the
% gradient at a vertex is computed by weighting the gradients of each
% triangle face sharing the vertex in proportion to the inverse distance
% from the center of each triangle to the vertex.
%
% [dFx,dFy] = TRIGRADIENT(TRI,X,Y,Z,'face') returns the numerical gradients
% that are constant over each triangle face. In this case, dFx and dFy are
% column vectors having M = size(TRI,1) elements.
%
% See also GRADIENT, DELAUNAY, TRIMESH, TRISURF, TRIPLOT.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-05-25, 2006-11-24

if nargin==4
   opt='vertex';
elseif nargin<4
	error('Not Enough Input Arguments.')
end
if ~ischar(opt)
   error('Fifth Input Argument Must be a String.')
end
x=x(:);	% convert input data into column vectors
y=y(:);
z=z(:);
xlen=length(x);
if ~isequal(xlen,length(y),length(z))
   error('X, Y, and Z Must Have the Same Number of Elements.')
end
if size(tri,2)~=3 || any(tri(:)<0) || any(tri(:)>xlen)
   error('TRI Must Be a Valid Triangulation of the Data in X, Y, Z.')
end

t=tri(:,[1 2 3 1]);
dy=diff(y(t),1,2);                           % [y2-y1 y3-y2 y1-y3]
dx=diff(x(t),1,2);                           % [x2-x1 x3-x2 x1-x3]
delta=-sum(dy(:,[2 3 1]).*x(tri),2);         % determinant
dxface=-sum(dy(:,[2 3 1]).*z(tri),2)./delta; % face gradient x
dyface= sum(dx(:,[2 3 1]).*z(tri),2)./delta; % face gradient y
if strncmpi(opt,'face',min(4,length(opt)))
   dfx=dxface;
   dfy=dyface;
   return
end
% vertex gradient requested
dfx=zeros(xlen,1);   % allocate space for results
dfy=zeros(xlen,1);
xc=sum(x(tri),2)/3;  % centroids of all triangles
yc=sum(y(tri),2)/3;

[vert,idx]=sort(tri(:));      % sort vertices in ascending order
tlen=size(tri,1);             % number of triangles
tn=[1:tlen 1:tlen 1:tlen]';   % triangle number for each vertex
tn=tn(idx);                   % shuffle triangle number to match vertices

last=find([diff(vert);1]);    % index of last vertex element in tn
first=1;

for k=1:xlen	% find gradient vertex by vertex

   r=tn(first:last(k));                   % triangles having k-th vertex
   first=last(k)+1;                       % beginning index for next vertex
   idist=1./hypot(x(k)-xc(r),y(k)-yc(r)); % inverse distances from vertex
   w=idist.'./sum(idist);                 % inverse distance weightings
   dfx(k)=w*dxface(r);
   dfy(k)=w*dyface(r);
end
   