
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Section


%mat file for saving mesh info
savename = 'mesh_WUS_Johnson2023';

%origin for Caresian coordinate system
origin = [34 -120];

%%
%load GPS data file
%format of columns:
%   lon       lat       Ve       Vn       Se       Sn      
data = load('Zeng_vels_augmented_PA_JDF.txt');


%% load creeping fault data

%If there are no creeping faults, set creeping_faults = [];

%format of columns:
% fault id, lon endpoint 1, lat endpoint 1, lon endpoint 2, lat endpoint 2, creep rate (mm/yr)
% NOTE: fault id is an integer, a unique identifier for each fault that section is considered to be a "continuous fault" 
creeping_faults = load('creeping_faults.txt');
%creeping_faults = [];
patchL = 15;  %maximum patch length


%% mesh domain size
%nominal node spacing
nom_node_spacing = 100;  %km
%longitude and latitude range of meshed domain
lon_range = [-127 -96];
lat_rage = [26 54];
%option to refine mesh in vicinity of GPS data
%specifiy level of refinement, integer from 0 to 4
%level 0 means no refinement, 4 is maximum refinement
refine_mesh = 0;


%% End Input Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath mesh2d
addpath tools




%convert from lon,lat to local Cartesian
xy_gps = (llh2local(data(:,1:2)',fliplr(origin)))';
Ve = data(:,3);  %east component of velocity
Vn = data(:,4);  %north component of velocity
Sige = data(:,5);  %standard deviation, east component
Sign = data(:,6);  %standard deviation, north component



%%
bxy = llh2local([lon_range;lat_rage],fliplr(origin))';
maxx = max(bxy(:,1)); 
minx = min(bxy(:,1));
maxy = max(bxy(:,2)); 
miny = min(bxy(:,2));

nx = round((maxx - minx)/nom_node_spacing);
ny = round((maxy-miny)/nom_node_spacing);
x=linspace(minx,maxx,nx);
y=linspace(miny,maxy,ny);
[X,Y]=meshgrid(x,y);
nodes = [X(:) Y(:)];

ind_gps = xy_gps(:,1)<maxx & xy_gps(:,1)>minx & xy_gps(:,2)<maxy & xy_gps(:,2)>miny;

%further refine in vicinity of GPS data by using nodes of voronoi cells in
%delaunay triangulation
[vx,vy]=voronoi(xy_gps(ind_gps,1),xy_gps(ind_gps,2));
vx = vx(:);
vy = vy(:);
[y,i]=unique(vx);
vx = vx(i);
vy = vy(i);

ind_vor = vx<maxx & vx>minx & vy<maxy & vy>miny;
nodes_refine =[vx(ind_vor) vy(ind_vor)];

if refine_mesh==4
    nodes = [nodes; xy_gps(ind_gps,:); nodes_refine];
elseif refine_mesh==3
    nodes = [nodes; xy_gps(ind_gps,:); nodes_refine(1:2:end,:)];
elseif refine_mesh==2
    nodes = [nodes; xy_gps(ind_gps,:); nodes_refine(1:3:end,:)];
elseif refine_mesh==1
    nodes = [nodes; xy_gps(ind_gps,:); nodes_refine(1:5:end,:)];
else
    nodes = [nodes; xy_gps(ind_gps,:)];
end

%get nodes on creeping faults
if ~isempty(creeping_faults)
    make_patches_creep
else
    node_creep = [];
    SegEnds = [];
    PatchEnds = [];
    PatchCreepRates = [];
end

cnt = size(nodes,1);
nodes = [nodes;node_creep];

tri = delaunay(nodes(:,1),nodes(:,2));

[edge] = tricon2(tri);
ebnd = edge(:,4) < +1;              %-- use bnd edge
conn = edge(ebnd,1:2);

if ~isempty(creeping_faults)
    edges = [conn;edge_creep+cnt];
else
   edges = conn;
end

[nodes,etri, tri,tnum] = smooth2(nodes,edges,tri);

%toss out GPS data outside of mesh domain
xy_gps(~ind_gps,:) = [];
Ve(~ind_gps) = [];
Vn(~ind_gps) = [];
Sige(~ind_gps) = [];
Sign(~ind_gps) = [];



%plot GPS data
figure;
patch('faces',tri(:,1:3),'vertices',nodes, ...
'facecolor','w', ...
'edgecolor',[.6,.6,.6]) ;
hold on; 
quiver(xy_gps(:,1),xy_gps(:,2),Ve,Vn,'r')
axis equal
if ~isempty(PatchEnds); plot(PatchEnds(:,[1 3])',PatchEnds(:,[2 4])','b','linewidth',2); end

patch_stuff=make_triangular_patch_stuff(tri,[nodes 0*nodes(:,1)]);

tri_areas = patch_stuff.area_faces;
tri_centroids = patch_stuff.centroids_faces;

save(savename, 'tri', 'nodes', 'tri_areas', 'tri_centroids','origin','PatchEnds','PatchCreepRates','xy_gps','Ve','Vn','Sige','Sign')    
    
   