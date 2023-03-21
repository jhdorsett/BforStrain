
%build Gcreep -- piecewise linear smoothing along faults
faultnums = unique(floor(Patch_id));
cnt = 0;
Gcreep = zeros(2*length(ind),size(pm,1));
GExx_creep = zeros(length(ind),size(pm,1));
GExy_creep = zeros(length(ind),size(pm,1));
GEyy_creep = zeros(length(ind),size(pm,1));

for k=1:length(faultnums)

   n = sum(floor(Patch_id)==faultnums(k));
   
   if n==1
   
        Gcreep(:,cnt+1) = [G1east_creep(:,cnt+1);G1north_creep(:,cnt+1)];  
        GExx_creep(:,cnt+1) = G1Exx_creep(:,cnt+1);  
        GExy_creep(:,cnt+1) = G1Exy_creep(:,cnt+1);  
        GEyy_creep(:,cnt+1) = G1Eyy_creep(:,cnt+1);  
       
   elseif n==2
       
       Gcreep(:,cnt+1) = [G1east_creep(:,cnt+1)-G2east_creep(:,cnt+1);G1north_creep(:,cnt+1)-G2north_creep(:,cnt+1)];  
       Gcreep(:,cnt+2) = [G2east_creep(:,cnt+1)+G1east_creep(:,cnt+2);G2north_creep(:,cnt+1)+G1north_creep(:,cnt+2)];  

       GExx_creep(:,cnt+1) = G1Exx_creep(:,cnt+1)-G2Exx_creep(:,cnt+1);  
       GExx_creep(:,cnt+2) = G2Exx_creep(:,cnt+1)+G1Exx_creep(:,cnt+2);  
       
       GExy_creep(:,cnt+1) = G1Exy_creep(:,cnt+1)-G2Exy_creep(:,cnt+1);  
       GExy_creep(:,cnt+2) = G2Exy_creep(:,cnt+1)+G1Exy_creep(:,cnt+2);  
       
       GEyy_creep(:,cnt+1) = G1Eyy_creep(:,cnt+1)-G2Eyy_creep(:,cnt+1);  
       GEyy_creep(:,cnt+2) = G2Eyy_creep(:,cnt+1)+G1Eyy_creep(:,cnt+2);  

   else
 
        Gcreep(:,cnt+1) = [G1east_creep(:,cnt+1)-G2east_creep(:,cnt+1);G1north_creep(:,cnt+1)-G2north_creep(:,cnt+1)];  
        GExx_creep(:,cnt+1) = G1Exx_creep(:,cnt+1)-G2Exx_creep(:,cnt+1);  
        GExy_creep(:,cnt+1) = G1Exy_creep(:,cnt+1)-G2Exy_creep(:,cnt+1);  
        GEyy_creep(:,cnt+1) = G1Eyy_creep(:,cnt+1)-G2Eyy_creep(:,cnt+1);  
 
       for j=cnt+2:cnt+(n-1)
            Gcreep(:,j) = [G2east_creep(:,j-1)+G1east_creep(:,j)-G2east_creep(:,j);G2north_creep(:,j-1)+G1north_creep(:,j)-G2north_creep(:,j)];  
            GExx_creep(:,j) = G2Exx_creep(:,j-1)+G1Exx_creep(:,j)-G2Exx_creep(:,j);  
            GExy_creep(:,j) = G2Exy_creep(:,j-1)+G1Exy_creep(:,j)-G2Exy_creep(:,j);  
            GEyy_creep(:,j) = G2Eyy_creep(:,j-1)+G1Eyy_creep(:,j)-G2Eyy_creep(:,j);  
            
       end     
       
        Gcreep(:,cnt+n) = [G2east_creep(:,cnt+n-1)+G1east_creep(:,cnt+n);G2north_creep(:,cnt+n-1)+G1north_creep(:,cnt+n)];  
        GExx_creep(:,cnt+n) = G2Exx_creep(:,cnt+n-1)+G1Exx_creep(:,cnt+n);  
        GExy_creep(:,cnt+n) = G2Exy_creep(:,cnt+n-1)+G1Exy_creep(:,cnt+n);  
        GEyy_creep(:,cnt+n) = G2Eyy_creep(:,cnt+n-1)+G1Eyy_creep(:,cnt+n);  
 
             
   end
   
   

   cnt = cnt+n;
   
end

Geast_creep = Gcreep(1:end/2,:);
Gnorth_creep = Gcreep(1+end/2:end,:);
