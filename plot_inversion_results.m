Exx_r = Exx_realizations;
Exy_r = Exy_realizations;
Eyy_r = Eyy_realizations;

%compute maximum shear strain and dilation rate for all realizations
max_shear_realizations = sqrt((Exx_r-Eyy_r).^2 + Exy_r.^2);
dilatation_realizations = Exx_r+Eyy_r;


%mean ane std of all realizations
Exx_mean = mean(Exx_r,2);
Exy_mean = mean(Exy_r,2);
Eyy_mean = mean(Eyy_r,2);

Exx_std = std(Exx_r,[],2);
Exy_std = std(Exy_r,[],2);
Eyy_std = std(Eyy_r,[],2);

mean_maxshear = sqrt((Exx_mean-Eyy_mean).^2 + Exy_mean.^2);
std_maxshear = std(max_shear_realizations,[],2);

mean_dilatation = Exx_mean + Eyy_mean;
std_dilatation = std(dilatation_realizations,[],2);



%% color plot of strain rates
figure
subplot(121)
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[min(mean_maxshear), max(mean_maxshear)]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData',mean_maxshear,...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(jet)
colorbar
title('Mean Maximum Shear Strain Rate (micro-strain/yr)')
set(gca,'ColorScale','log')
caxis([0.01 1])

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end


subplot(122)
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[min(2*std_maxshear), max(2*std_maxshear)]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData',2*std_maxshear,...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(jet)
colorbar
title('Two standard deviations (micro-strain/yr)')
set(gca,'ColorScale','log')
caxis([0.01 1])

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

drawnow


%% dilatation rates
figure
subplot(121)
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[min(mean_dilatation), max(mean_dilatation)]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData',mean_dilatation,...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(jet)
colorbar
title('Mean Dilatation Rate (micro-strain/yr)')
caxis([-.1 .1])

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end


subplot(122)
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[min(2*std_dilatation), max(2*std_dilatation)]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData',2*std_dilatation,...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
load cmap
colormap(flipud(cmap))
colorbar
title('Two standard deviations (micro-strain/yr)')
caxis([-.1 .1])

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

drawnow

%% plot strain rates exceeding 2-sigma uncertainties
ind = mean_maxshear>2*std_maxshear;

figure
subplot(121)
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title(' Maximum Shear Strain Rate exceeding 2\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end


ind = mean_maxshear>1*std_maxshear;
subplot(122)

hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title(' Maximum Shear Strain Rate exceeding 1\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

drawnow


%% plot ratio of std to mean
rat = abs(mean_maxshear)./std_maxshear;

figure
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', rat,...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
caxis([.5 5])
colorbar
title(' STD divided by Mean Max Shear')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

drawnow

%% plot ratio of std to mean
rat = abs(mean_dilatation)./std_dilatation;

figure
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', rat,...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
caxis([.5 5])
colorbar
title(' STD divided by Mean Dilatation')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

drawnow



%% plot dilatancy exceeding 2-sigma uncertainties
ind = abs(mean_dilatation)>2*std_dilatation;

figure
subplot(121)
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title(' Dilatation Rate exceeding 2\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end



ind = abs(mean_dilatation)>1*std_dilatation;
subplot(122)

hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title(' Dilatation Rate exceeding 1\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end


drawnow




%plot strain rate component exceeding 2-sigma uncertainties
ind = abs(Exx_mean)>2*Exx_std |  abs(Exy_mean)>2*Exy_std |  abs(Eyy_mean)>2*Eyy_std  ;

figure
subplot(121)
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title(' One Strain Rate Component exceeding 2\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end



ind = abs(Exx_mean)>1*Exx_std |  abs(Exy_mean)>1*Exy_std |  abs(Eyy_mean)>1*Eyy_std  ;
subplot(122)

hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title(' One Strain Rate Component exceeding 1\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end



drawnow



%% principal directions
for j=1:size(Exx_r,2)

        for k=1:size(Exx_r(:,j),1)

            E = [Exx_r(k,j) Exy_r(k,j); Exy_r(k,j) Eyy_r(k,j)];
            [vec,val] = eig(E);

            minVecs(:,k,j) = vec(:,1);
            maxVecs(:,k,j) = vec(:,2);
            minvals(k,j) = val(1,1);
            maxvals(k,j) = val(2,2);

        end
        
end




%show where principal strain rates exceed zero
ind1  = (abs(mean(minvals,2)) > std(minvals,[],2)) | (abs(mean(maxvals,2)) > std(maxvals,[],2));
ind2  = (abs(mean(minvals,2)) > 1.5*std(minvals,[],2)) | (abs(mean(maxvals,2)) > 1.5*std(maxvals,[],2));
ind3  = (abs(mean(minvals,2)) > 2*std(minvals,[],2)) | (abs(mean(maxvals,2)) > 2*std(maxvals,[],2));


figure
subplot(131)
hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind3),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title('One Principal strain Rate exceeding 2\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end



subplot(132)

hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind2),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title(' One Principal Strain Rate exceeding 1.5\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end



subplot(133)

hh = trisurf(tri,nodes(:,1),nodes(:,2),0*nodes(:,2));
set(gca,'CLim',[0, 2]);
set(hh,'FaceColor','flat',...
       'FaceVertexCData', double(ind1),...
       'CDataMapping','scaled');
set(hh,'edgecolor','none');
axis equal
view(2)
colormap(flipud(hot))

title(' One Principal Strain Rate exceeding 1\sigma uncertainties')

hold on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

drawnow

vecscale = 5;
figure
hold on
skip=2;
% for k=1:skip:size(Exx_r,1)
% 
%            angles_max = squeeze(atan(maxVecs(2,k,:)./maxVecs(1,k,:))*180/pi);
%             [dip, p_value, xlow, xup] = HartigansDipSignifTest(angles_max, 50);
%                     if p_value<0.5 & mean(abs(angles_max))>45
%                         angles_max(angles_max<0) = 180 + angles_max(angles_max<0);
%                     end
%             angles_min = squeeze(atan(minVecs(2,k,:)./minVecs(1,k,:))*180/pi);
%             [dip, p_value, xlow, xup] = HartigansDipSignifTest(angles_min, 50);
%                     if p_value<0.5  & mean(abs(angles_min))>45
%                         angles_min(angles_min<0) = 180 + angles_min(angles_min<0);
%                     end
%                     
%             if abs(mean(minvals(k,:))) > abs(mean(maxvals(k,:)))
%                  %plot maximum shortening rate direction
%                     angles = angles_min;
%                     sign_val = sign(mean(minvals(k,:)));
%                     mean_angle2 = mean(angles_max); 
%                     sign_val2 = sign(mean(maxvals(k,:)));
% 
%                 
%             else
%                 angles = angles_max;      
%                 sign_val = sign(mean(maxvals(k,:)));
%                 mean_angle2 = mean(angles_min); 
%                 sign_val2 = sign(mean(minvals(k,:)));
%             end
%             
%               std_angles = std(angles);
% 
%               
%             if std_angles<40
%             
%           
%             
%              %if abs(mean(minvals(k,:)))>2*std(minvals(k,:)) | abs(mean(maxvals(k,:)))>2*std(maxvals(k,:))   
%                 a1 = pi/180*(mean(angles) + 2*std_angles);
%                 a2 = pi/180*(mean(angles) - 2*std_angles);
%                 t = linspace(a1,a2,128);
% 
%                 x0 = [0 vecscale*cos(t) 0]+tri_centroids(k,1);
%                 z0 = [0 vecscale*sin(t) 0]+tri_centroids(k,2);
%                 if sign_val<0
%                     patch( x0, z0, 'r','FaceAlpha',.3 ,'EdgeColor','none');
%                 else
%                     patch( x0, z0, 'b' ,'FaceAlpha',.3,'EdgeColor','none');
%                 end
% 
%                 x0 = [0 -vecscale*cos(t) 0]+tri_centroids(k,1);
%                 z0 = [0 -vecscale*sin(t) 0]+tri_centroids(k,2);
% 
%                 if sign_val<0
%                     patch( x0, z0, 'r', 'FaceAlpha',.3,'EdgeColor','none');
%                 else
%                     patch( x0, z0, 'b', 'FaceAlpha',.3,'EdgeColor','none');
%                  end
% 
%             end
% 
%         
%              
%             %plot tick in other direction
%             xs = [tri_centroids(k,1)+vecscale*cos(mean_angle2*pi/180) tri_centroids(k,1)-vecscale*cos(mean_angle2*pi/180)];
%             ys = [tri_centroids(k,2)+vecscale*sin(mean_angle2*pi/180) tri_centroids(k,2)-vecscale*sin(mean_angle2*pi/180)];
% 
%     if sign_val2>0
%             plot(xs,ys,'b')
%     else
%             plot(xs,ys,'r')
%     end
%    
% end
mean_minvals = mean(minvals,2);
mean_maxvals = mean(maxvals,2);

for k=1:skip:size(Exx_r,1)

    %plot the larger magnitude principal direction
    
    if abs(mean_minvals(k))>abs(mean_maxvals(k))
        
        bigvals = mean_minvals(k);
        smallvals = mean_maxvals(k);
        
        bigVecs = squeeze(minVecs(:,k,:));
        smallVecs = squeeze(maxVecs(:,k,:));
        

                
    else
        
        bigvals = mean_maxvals(k);
        smallvals = mean_minvals(k);        
    
        bigVecs = squeeze(maxVecs(:,k,:));
        smallVecs = squeeze(minVecs(:,k,:));

    end
    
    mean_smallVecs = mean(smallVecs,2);
    angles = squeeze(atan(bigVecs(2,:)./bigVecs(1,:))*180/pi);
    
    
    %determin if bimodal. If so, add 180 to negative values
    [dip, p_value, xlow, xup] = HartigansDipSignifTest(angles, 50);
    if p_value<0.5 & mean(abs(angles))>45
        angles(angles<0) = 180 + angles(angles<0);
    end
    
    std_angles = std(angles);
    
    if std_angles<40
        
        a1 = pi/180*(mean(angles) + 2*std_angles);
        a2 = pi/180*(mean(angles) - 2*std_angles);
        t = linspace(a1,a2,128);



        x0 = [0 vecscale*cos(t) 0]+tri_centroids(k,1);
        z0 = [0 vecscale*sin(t) 0]+tri_centroids(k,2);    
        if bigvals<0
             patch( x0, z0, 'r','FaceAlpha',.5 ,'EdgeColor','none');
        else
             patch( x0, z0, 'b','FaceAlpha',.5 ,'EdgeColor','none');
        end    

        x0 = [0 -vecscale*cos(t) 0]+tri_centroids(k,1);
        z0 = [0 -vecscale*sin(t) 0]+tri_centroids(k,2);

        if bigvals<0
           patch( x0, z0, 'r','FaceAlpha',.5 ,'EdgeColor','none');
        else
            patch( x0, z0, 'b','FaceAlpha',.5 ,'EdgeColor','none');
        end



        %plot smaller principal direction as a  line
        vx = [tri_centroids(k,1)-mean_smallVecs(1)*vecscale*abs(smallvals/bigvals) tri_centroids(k,1)+mean_smallVecs(1)*vecscale*abs(smallvals/bigvals)];
        vy = [tri_centroids(k,2)-mean_smallVecs(2)*vecscale*abs(smallvals/bigvals) tri_centroids(k,2)+mean_smallVecs(2)*vecscale*abs(smallvals/bigvals)];


        if smallvals<0
            plot(vx,vy,'r')
        else
            plot(vx,vy,'b')
        end

    end

    style(k) = (bigvals+smallvals)/(abs(bigvals)+abs(smallvals));
    
end
grid on

if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

axis equal

title('2\sigma error wedges, maximum shortening rate direction')



drawnow


%% plot only tick marks


vecscale = 5;
figure
hold on

for k=1:2:size(Exx_r,1)

    
        angles = squeeze(atan(minVecs(2,k,:)./minVecs(1,k,:))*180/pi);
        
         %determin if bimodal. If so, add 180 to negative values
            [dip, p_value, xlow, xup] = HartigansDipSignifTest(angles, 50);
            if p_value<0.5 & mean(abs(angles))>45
                angles(angles<0) = 180 + angles(angles<0);
            end
            
            
        angle = mean(angles);
 
            
    xs = [tri_centroids(k,1)+vecscale*cos(angle*pi/180) tri_centroids(k,1)-vecscale*cos(angle*pi/180)];
    ys = [tri_centroids(k,2)+vecscale*sin(angle*pi/180) tri_centroids(k,2)-vecscale*sin(angle*pi/180)];

    
    plot(xs,ys,'b')


        
end
grid on


if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

axis equal
title('Maximum shortening rate direction')




%plot style

figure
hold on
scatter(tri_centroids(1:skip:end,1),tri_centroids(1:skip:end,2),50,-style(1:skip:end),'fill')
grid on
if exist('SegEnds'); if ~isempty(SegEnds); plot(SegEnds(:,[1 3])',SegEnds(:,[2 4])','k'); end; end

axis equal
colormap(jet)
colorbar
title('Strain Rate Style')




