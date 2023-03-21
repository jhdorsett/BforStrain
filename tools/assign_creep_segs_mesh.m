
creep_rate = creeping_faults(:,5);

make_patches_creep

[lia1,locb1] = ismember(pm_seg_num,1:length(creep_rate));
pm_creep_rate = creep_rate(locb1);

pm_id = GeoBounds(locb1,1);

pm = pm(pm_creep_rate>0,:);
pm_id = pm_id(pm_creep_rate>0);
pm_creep_rate = pm_creep_rate(pm_creep_rate>0);


%make Garlock left-lateral (negative)
ind = [75:76 4:15 91:98];
pm_creep_rate(ind) = -pm_creep_rate(ind);


segends1 = [pm(:,6)+ pm(:,1)/2.*cos((90-pm(:,5))*pi/180) pm(:,7)+ pm(:,1)/2.*sin((90-pm(:,5))*pi/180)];
segends2 = [pm(:,6)- pm(:,1)/2.*cos((90-pm(:,5))*pi/180) pm(:,7)- pm(:,1)/2.*sin((90-pm(:,5))*pi/180)];

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

%edge = edge(1:350,:);

%node_creep = node_creep(1:max(max(edge_creep)),:);
