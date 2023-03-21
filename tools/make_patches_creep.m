

%convert lat,long to x,y
x=creeping_faults(:,2:5);
llhx = x(:,1:2);
x1 = llh2local(llhx', fliplr(origin))';
llhx = x(:,3:4);
x2 = llh2local(llhx', fliplr(origin))';
SegEnds=[x1 x2];

rempatch = SegEnds(:,1)<minx | SegEnds(:,1)>maxx | SegEnds(:,2)<miny | SegEnds(:,2)>maxy;
rempatch = rempatch | SegEnds(:,3)<minx | SegEnds(:,3)>maxx | SegEnds(:,4)<miny | SegEnds(:,4)>maxy;
SegEnds(rempatch,:) = [];
creeping_faults(rempatch,:) = []; 

PatchEnds = [];
PatchCreepRates = [];
Patch_id = [];
for k=1:size(SegEnds,1)

    patchlength = sqrt((SegEnds(k,3)-SegEnds(k,1)).^2 + (SegEnds(k,4)-SegEnds(k,2)).^2);
    numpatch = ceil(patchlength/patchL);
    
    xs = linspace(SegEnds(k,1),SegEnds(k,3),numpatch+1);
    ys = linspace(SegEnds(k,2),SegEnds(k,4),numpatch+1);
    
    PatchEnds = [PatchEnds; [xs(1:end-1)' ys(1:end-1)' xs(2:end)' ys(2:end)']];
    PatchCreepRates = [ PatchCreepRates; creeping_faults(k,6)*ones(numpatch,1)];
    Patch_id = [Patch_id; creeping_faults(k,6)*ones(numpatch,1)];
    
end


segends1 = PatchEnds(:,1:2);
segends2 = PatchEnds(:,3:4);
node_creep = zeros(size(segends1,1)+size(segends2,1),2);

node_creep(1:2:end) = segends1;
node_creep(2:2:end) = segends2;

edge_creep = [(1:2:size(node_creep,1))' (2:2:size(node_creep,1))'];



%need to remove repeated nodes because
%meshing algorithm doesn't like repeats
node_creep = round(node_creep,4);
[C,ia,ic] = unique(node_creep,'rows','stable');
node_creep = C;


%need to renumber nodes values in edges to be consistent 
%with unique node values
new_edge = edge_creep;
for k=1:length(ic)
    new_edge(edge_creep==k)=ic(k);
end

edge_creep = new_edge;





 