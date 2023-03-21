clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

%load the mat file generated in setup_mesh
load mesh_wus_creep

%filename for saving Greens functions in mat file
saveGreens = 'test_Greens'; 

%% No need to modify below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%assign data to triangles
find_triangles_data

%take median of data in a triangle and assign this value to the centroid
assign_data_centroids


%creep Greens functions (piece-wise linear slip gradients)
refine = 1;  %creep discretization (km)

if ~isempty(PatchEnds)
    lengths = sqrt((PatchEnds(:,3)-PatchEnds(:,1)).^2 + (PatchEnds(:,4)-PatchEnds(:,2)).^2);
    angle = atan2(PatchEnds(:,4)-PatchEnds(:,2),PatchEnds(:,3)-PatchEnds(:,1));
    strike=90-angle*180/pi;
    centers=[(PatchEnds(:,1)+PatchEnds(:,3))./2 (PatchEnds(:,2)+PatchEnds(:,4))./2];
    pm = [lengths 10^6*ones(length(lengths),1) 10^6*ones(length(lengths),1) 90*ones(length(lengths),1) strike centers(:,1) centers(:,2)];
    [G1east_creep,G2east_creep,G1north_creep,G2north_creep,...
        G1Exx_creep,G2Exx_creep,G1Exy_creep,G2Exy_creep,...
        G1Eyy_creep,G2Eyy_creep] = make_dispG_piecewise(pm,tri_centroids(:,1:2),refine);
end

nu = .25;  %Poisson's ratio
[Ge_x,Ge_y,Gn_x,Gn_y,GExx_x,GExy_x,GEyy_x,GExx_y,GExy_y,GEyy_y] = buildG_PointForce(tri_centroids,nodes,nu);

fprintf('\n \n %s \n \n', 'Saving Greens Functions')

save(saveGreens)
