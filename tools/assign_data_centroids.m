Ve_centroids = nan(size(tri,1),1);
Vn_centroids = nan(size(tri,1),1);

Sige_centroids = nan(size(tri,1),1);
Sign_centroids = nan(size(tri,1),1);


for k=1:size(tri,1) 
    ind = (data_tri == k);
    if sum(ind)>0
        
        %weighted mean
        W = 1./Sige(ind);
        Ve_centroids(k) = sum(W.*Ve(ind)/sum(W));
        
        W = 1./Sign(ind);
        Vn_centroids(k) = sum(W.*Vn(ind)/sum(W));
       

        %propagate errors
        Sige_centroids(k) = sqrt(sum(Sige(ind).^2));
        Sign_centroids(k) = sqrt(sum(Sign(ind).^2));

               
    end
end

V_centroids = sqrt(Ve_centroids.^2 + Vn_centroids.^2);

figure;
patch('faces',tri(:,1:3),'vertices',nodes, ...
'facecolor','w', ...
'edgecolor',[.6,.6,.6]) ;
hold on;

for k=1:size(tri,1)
    if num_data_tri(k)~=0
        fill(nodes(tri(k,:),1),nodes(tri(k,:),2),V_centroids(k))
    end
end
axis equal
colorbar
colormap(jet)

title('observed velocities in each triangle')
set(gca,'fontsize',15)
