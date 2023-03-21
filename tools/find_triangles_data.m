
data_tri = zeros(size(xy_gps,1),1);
num_data_tri = zeros(size(tri,1),1);
for k=1:size(tri,1)
    
    IN = inpolygon(xy_gps(:,1),xy_gps(:,2),nodes(tri(k,:),1),nodes(tri(k,:),2));

    data_tri(IN) = k;
    num_data_tri(k) = sum(IN);

end

figure;
patch('faces',tri(:,1:3),'vertices',nodes, ...
'facecolor','w', ...
'edgecolor',[.6,.6,.6]) ;
hold on;
colorbar
colormap(jet)

for k=1:size(tri,1)
    if num_data_tri(k)~=0
        fill(nodes(tri(k,:),1),nodes(tri(k,:),2),num_data_tri(k))
    end
end
title('number of data points in each triangle')
set(gca,'fontsize',15)

axis equal



