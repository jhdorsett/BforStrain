function [G1east,G2east,G1north,G2north,G1Exx,G2Exx,G1Exy,G2Exy,G1Eyy,G2Eyy] = make_dispG_piecewise(pm,xystats,refine)


xloc=[xystats';zeros(1,size(xystats',2))];

npatches=size(pm,1);

G1east = zeros(size(xystats,1),npatches);
G2east = zeros(size(xystats,1),npatches);
G1north = zeros(size(xystats,1),npatches);
G2north = zeros(size(xystats,1),npatches);
G1Exx = zeros(size(xystats,1),npatches);
G2Exx = zeros(size(xystats,1),npatches);
G1Exy = zeros(size(xystats,1),npatches);
G2Exy = zeros(size(xystats,1),npatches);
G1Eyy = zeros(size(xystats,1),npatches);
G2Eyy = zeros(size(xystats,1),npatches);


for k=1:npatches
    
    %divide into small segments 
    nhe = ceil(pm(k,1)/refine);
    pf = patchfault(pm(k,:),nhe,1);
    
    east1 = zeros(size(xystats,1),1);
    east2 = zeros(size(xystats,1),1);
    north1 = zeros(size(xystats,1),1);
    north2 = zeros(size(xystats,1),1);
    Exx1 = zeros(size(xystats,1),1);
    Exx2 = zeros(size(xystats,1),1);
    Exy1 = zeros(size(xystats,1),1);
    Exy2 = zeros(size(xystats,1),1);
    Eyy1 = zeros(size(xystats,1),1);
    Eyy2 = zeros(size(xystats,1),1);
    
    for j=1:nhe
        
        
        m1=[pf(j,:) -1 0 0]';
   
        [U1,D,S]=disloc3d(m1,xloc,1,.25);
       
        east1= east1 + U1(1,:)';
        north1 = north1 + U1(2,:)';

        east2= east2 + (j-1)/nhe*U1(1,:)';
        north2 = north2 + (j-1)/nhe*U1(2,:)';

        
        Exx = D(1,:)';
        Exy = .5*(D(2,:)+D(4,:))';
        Eyy = D(5,:)';

        Exx1 = Exx1 + Exx;
        Exy1 = Exy1 + Exy;
        Eyy1 = Eyy1 + Eyy;
        

        Exx2 = Exx2 + (j-1)/nhe*Exx;
        Exy2 = Exy2 + (j-1)/nhe*Exy;
        Eyy2 = Eyy2 + (j-1)/nhe*Eyy;

        
    end
    
    G1east(:,k) = east1;
    G2east(:,k) = east2;

    G1north(:,k) = north1;
    G2north(:,k) = north2;

    G1Exx(:,k) = Exx1;
    G1Exy(:,k) = Exy1; 
    G1Eyy(:,k) = Eyy1;

    G2Exx(:,k) = Exx2;
    G2Exy(:,k) = Exy2; 
    G2Eyy(:,k) = Eyy2;

    
    disp(['Completed ' num2str(k/npatches*100) '% of creeping fault Greens function calculations']) 
end

