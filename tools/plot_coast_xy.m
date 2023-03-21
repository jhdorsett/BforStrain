function plot_coast_xy(bbox,origin)
%bbox [minlon minlat;maxlon maxlat]

S = shaperead('ne_10m_coastline','BoundingBox',bbox);



for k=1:length(S)

    xy = llh2local([S(k).X; S(k).Y],[origin(1) ;origin(2)])';
    
    
    plot(xy(:,1),xy(:,2),'k')
    
end
