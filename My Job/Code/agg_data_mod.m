function [Agg,Ar,msk,Area]=agg_data_mod(data,basins,id)
% ===========================================================
%      function Agg=agg_data_mod(data,basins)
%
%  This function returns the aggregated time series from
%  grids based data (this version is meant  for 0.5 x 0.5 degree)

%  input:
%        data:    in cell format and year infirst column, month in second
%                 3rd column should conatin the filed (360 x 720)
%        basins:  an extracted shapefile in matlab using shaperead function
%        ids:      the row number belonging to the catchment
%   output:
%        Agg : the aggregated time series
%        Ar : area of the basin

% Mohammad J. Tourian
% tourian@gis.uni-stuttgart.de
% January 2012
% ===========================================================

[s1 s2]=size(data);
 Lat_ref=(89.75:-.5:-89.75)'*ones(1,720);
 Lon_ref=((-179.75:.5:179.75)'*ones(1,360))';
 
%  Lat_ref=(88.75:-2.5:-88.75)'*ones(1,144);
% Lon_ref=((-178.75:2.5:178.75)'*ones(1,72))';



R=6378137;
for i=1:length(id)
    ids=id(i);
    
    
    [m n]=find(Lon_ref>=basins(ids,1).BoundingBox(1,1)&Lon_ref<=basins(ids,1).BoundingBox(2,1)&Lat_ref>=basins(ids,1).BoundingBox(1,2)&Lat_ref<=basins(ids,1).BoundingBox(2,2));
    p=find(Lon_ref>=basins(ids,1).BoundingBox(1,1)&Lon_ref<=basins(ids,1).BoundingBox(2,1)&Lat_ref>=basins(ids,1).BoundingBox(1,2)&Lat_ref<=basins(ids,1).BoundingBox(2,2));
%     lat=basins(ids,1).Lat(find(isnan(basins(ids,1).Lat)==0));
%     lon=basins(ids,1).Lon(find(isnan(basins(ids,1).Lon)==0));

        lat=basins(ids,1).Lat;
    lon=basins(ids,1).Lon;
    
    INbasin=inpolygon(Lon_ref,Lat_ref,lon,lat);
    
    INp=find(INbasin==1);
    Otp=find(INbasin==0);
    
    Ar_b=max(areaint(lat,lon)*4*pi*R^2);
    box(1,:)=Lat_ref(p)-.25;
    box(2,:)=Lat_ref(p)+.25;
    box(3,:)=Lon_ref(p)-.25;
    box(4,:)=Lon_ref(p)+.25;
    
    
    [b1,b2]=size(box);
    for ii=1:b2
        xbox(ii,:)=[box(3,ii) box(3,ii) box(4,ii) box(4,ii) box(3,ii) ];
        ybox(ii,:)=[box(1,ii) box(2,ii) box(2,ii) box(1,ii) box(1,ii)];
    end
    
    Area=nan(360,720);
    
    for j=1:b2
        IN = inpolygon(xbox(j,:)',ybox(j,:)',lon,lat);
        f=find(IN==0);
        g=find(IN==1);
        if length(f)~=5
            if length(g)==5
                Area(m(j),n(j))=areaint(ybox(j,:),xbox(j,:))*4*pi*R^2;
            else
                [xi, yi] = polyxpoly(lon, lat, xbox(j,:), ybox(j,:));
                latemp=[ybox(j,g) yi'];
                lonemp=[xbox(j,g) xi'];
                Area(p(j))=areaint(latemp,lonemp)*4*pi*R^2;
                
                clear latemp lonemp
            end
        else
            
            Area(m(j),n(j))=areaint(ybox(j,:),xbox(j,:))*4*pi*R^2;
        end
        
    end
    Ar(i)=max(areaint(lat,lon)*4*pi*R^2);
    msk=Area;
    msk(Otp)=NaN;
    
          for j=1:s1
            mult=data{j,3}(INp).*Area(INp);
            Agg(j,i+2)=nansum(mult(:))/Ar(i);
            Agg(j,1)=data{j,1};Agg(j,2)=data{j,2};
        end
    
    clear m n box
end


