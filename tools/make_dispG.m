function [Geast,Gnorth,GExx,GExy,GEyy] = make_dispG(pm,xystats)

xloc=[xystats';zeros(1,size(xystats',2))];

npatches=size(pm,1);

Geast = zeros(size(xystats,1),npatches);
Gnorth = zeros(size(xystats,1),npatches);
GExx = zeros(size(xystats,1),npatches);
GExy = zeros(size(xystats,1),npatches);
GEyy = zeros(size(xystats,1),npatches);


for k=1:npatches
   m1=[pm(k,:) -1 0 0]';
   
   [U1,D,S,flag]=disloc3d(m1,xloc,1,.25);
   
   Geast(:,k)=U1(1,:)';
   Gnorth(:,k)=U1(2,:)';
   
   
    Exx = D(1,:)';
    Exy = .5*(D(2,:)+D(4,:))';
    Eyy = D(5,:)';

    GExx(:,k) = Exx;
    GExy(:,k) = Exy; 
    GEyy(:,k) = Eyy;

   
end

