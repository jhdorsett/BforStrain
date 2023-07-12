function [Ge_x, Ge_y, Gn_x, Gn_y,GExx_x,GExy_x,GEyy_x,GExx_y,GExy_y,GEyy_y,Omega_x,Omega_y] = buildG_PointForce(xy_gps,nodes,nu)

Ge_x = zeros(size(xy_gps,1),size(nodes,1));
Ge_y = zeros(size(xy_gps,1),size(nodes,1));
Gn_x = zeros(size(xy_gps,1),size(nodes,1));
Gn_y = zeros(size(xy_gps,1),size(nodes,1));
GExx_x = zeros(size(xy_gps,1),size(nodes,1));
GEyy_x = zeros(size(xy_gps,1),size(nodes,1));
GExy_x = zeros(size(xy_gps,1),size(nodes,1));
GExx_y = zeros(size(xy_gps,1),size(nodes,1));
GEyy_y = zeros(size(xy_gps,1),size(nodes,1));
GExy_y = zeros(size(xy_gps,1),size(nodes,1));
Omega_x = zeros(size(xy_gps,1),size(nodes,1));
Omega_y = zeros(size(xy_gps,1),size(nodes,1));


for k=1:size(nodes,1)
    
    [Ge_x(:,k), Ge_y(:,k), Gn_x(:,k), Gn_y(:,k), GExx_x(:,k),GExy_x(:,k),...
        GEyy_x(:,k),GExx_y(:,k),GExy_y(:,k),GEyy_y(:,k),Omega_x(:,k),Omega_y(:,k)] = PointForce(nodes(k,1),nodes(k,2),xy_gps(:,1),xy_gps(:,2),nu);
    
     disp(['Completed ' num2str(k/size(nodes,1)*100) '% of body force Greens function calculations']) 
     
end
