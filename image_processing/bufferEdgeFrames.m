function [imageOut] = bufferEdgeFrames(imageIn,imageThreshold)

s=size(imageIn);
center=s./2;

% Get region properties of largest blob
props = regionprops(imageIn>imageThreshold,'Centroid','Area');

if ~isempty(props)
    [v,i]=max([props.Area]);
    props=props(i);
    cen=round(props.Centroid);
    xDev = cen(1) - center(2);
    yDev = cen(2) - center(1);
    
    if abs(xDev) > 0.1*s(2) || abs(yDev) > 0.1*s(1)
        tmpIm = uint8(zeros(s.*2));
        tmp_center = round(size(tmpIm)./2);
        
        tmpIm(tmp_center(1)-center(1)-yDev+1:tmp_center(1)+center(1)-yDev,...
            tmp_center(2)-center(2)+1-xDev:tmp_center(2)+center(2)-xDev) = imageIn;
        
        imageOut = tmpIm(tmp_center(1)-center(1)+1:tmp_center(1)+center(1),...
            tmp_center(2)-center(2)+1:tmp_center(2)+center(2));
    else
        imageOut = imageIn;
    end
else
    imageOut = imageIn;
end
    


