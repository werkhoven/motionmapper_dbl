function [Xs,Ys,angles,areas,parameters,framesToCheck,svdskipped,areanorm] = ...
    alignImages_Radon_parallel_avi(file_path,startImage,finalImage,image_path,parameters)
%alignImages_Radon_parallel_avi runs the alignment and segmentation routines on a .avi file
%   and saves the output files to a directorty (called by ../runAlignment.m)
%
%   Input variables:
%
%       file_path -> avi file to be analyzed
%       startImage -> first frame of the avi file to be analyzed
%       finalImage -> last frame of the avi file to be analyzed
%       image_path -> path to which files are saved
%       parameters -> struct containing parameters
%
%
%   Output variables:
%
%       Xs -> alignment x translations
%       Ys -> alignment y translations
%       angles -> alignment rotations
%       areas -> segmented areas after segmentation
%       framesToCheck -> frames where a large rotation occurs in a single 
%                           frame.  This might signal a 180 degree rotation
%                           error
%       svdskipped -> blank frames where alignment is skipped
%       areanorm -> normalizing factor for fly size
%
% (C) Gordon J. Berman, 2014
%     Princeton University

%% Set processing parameters

    warning off MATLAB:polyfit:RepeatedPointsOrRescale;
    warning off MATLAB:audiovideo:aviinfo:FunctionToBeRemoved;
    
    readout = 100;
    nDigits = 8;
    
    spacing = parameters.alignment_angle_spacing;
    pixelTol = parameters.pixelTol;
    minArea = parameters.minArea;
    asymThreshold = parameters.asymThreshold;
    symLine = parameters.symLine;
    initialPhi = parameters.initialPhi;
    dilateSize = parameters.dilateSize;
    cannyParameter = parameters.cannyParameter;
    imageThreshold = parameters.imageThreshold;
    maxAreaDifference = parameters.maxAreaDifference;
    segmentationOff = parameters.segmentationOff;
    basisImage = parameters.basisImage;
    bodyThreshold = parameters.bodyThreshold;
    numProcessors = parameters.numProcessors;
    rangeExtension = parameters.rangeExtension;
    
    %Choose starting and finishing images
    
    vidObj = VideoReader(file_path);
    [~,fLabel,~] = fileparts(file_path);
    nFrames = vidObj.NumberOfFrames;
    
    if isempty(startImage)
        startImage = 1000;
    end
    
    if isempty(finalImage)
        if nFrames < 361000
            finalImage = nFrames;
        else
            finalImage = 361000;
        end
    end
    
        
    segmentationOptions.imageThreshold = imageThreshold;
    segmentationOptions.cannyParameter = cannyParameter;
    segmentationOptions.dilateSize = dilateSize;
    segmentationOptions.minArea = minArea;
    segmentationOptions.spacing = spacing;
    segmentationOptions.pixelTol = pixelTol;
    segmentationOptions.maxAreaDifference = maxAreaDifference;
    segmentationOptions.segmentationOff = segmentationOff;
    segmentationOptions.asymThreshold = asymThreshold;
    segmentationOptions.symLine = symLine;
    
%% Load centroid data


    camRes = [1024 1280];
    refpad = 100;
    [fDir,fLabel,~] = fileparts(file_path);
    camNum = str2double(fLabel(strfind(fLabel,'cam')+3));
    [cenPath] = getHiddenMatDir(fDir,'ext','.txt','Keyword','CentroidData');
    cenDat = dlmread(cenPath{:});
    
    tmp = ones(length(cenDat)/12,2);
    tmp(:,1) = cenDat(mod(1:length(cenDat),12)==camNum*3-1);
    if camNum<4
        tmp(:,2) = cenDat(mod(1:length(cenDat),12)==camNum*3);
    else
        tmp(:,2) = cenDat(mod(1:length(cenDat),12)==0);
    end
    cenDat = tmp;
    clearvars tmp
 
%% Load reference image or build one from scratch

refPath = [fDir '\' fLabel(1:end-5) '_ref' num2str(camNum) '.png'];
if exist(refPath,'file')
    
        ref = imread(refPath);
        ref2 = uint8(zeros(camRes+refpad*2));
        ref2(refpad+1:refpad+camRes(1),refpad+1:refpad+camRes(2)) = ref;
        ref = ref2;
        clear ref2
        
else

    % define grid spacing for reference building
    ref = NaN(camRes);
    stp_sz = 4;
    xStp = length(stp_sz:stp_sz:camRes(2));
    yStp = length(stp_sz:stp_sz:camRes(1));
    iFr = NaN(xStp*yStp,3);
    ct = 0;
    
    % Find frames matching the desired spacing
    hwb = waitbar(0,'finding target frames for reference tiling');
    for i = 1:xStp
        for j = 1:yStp
            
            ct = ct+1;
            targetX = cenDat(:,1) == i*stp_sz;
            targetY = cenDat(:,2) == j*stp_sz;
            onTarget = find(targetX & targetY,1);
            
            if ~isempty(onTarget)
                iFr(ct,1) = onTarget;
                iFr(ct,2) = i*stp_sz;
                iFr(ct,3) = j*stp_sz;
            else
                dx = abs(cenDat(:,1) - i*stp_sz)<3;
                dy = abs(cenDat(:,2) - j*stp_sz)<3;
                nearTarget = find(dx & dy,1);
                if ~isempty(nearTarget) && ~any(iFr(:,2)==cenDat(nearTarget,1)&iFr(:,3)==cenDat(nearTarget,2))
                    iFr(ct,1) = nearTarget;
                    iFr(ct,2) = cenDat(nearTarget,1);
                    iFr(ct,3) = cenDat(nearTarget,2);
                end               
            end
            
            hwb = waitbar(ct/(xStp*yStp));
            
        end
    end
    
    delete(hwb);
    iFr(isnan(iFr(:,1)),:)=[];
    hwb = waitbar(0,'','Name','sampling images to build reference');
    
    % Build reference image from queried frames
    for i = 1:length(iFr)
        
        hwb = waitbar(i/length(iFr),hwb,['frame ' num2str(i) ' of ' num2str(length(iFr))]);
        
        % grab current frame
        if iFr(i,1) <= nFrames
        tmp_frame = read(vidObj,iFr(i,1));
        tmp_frame(45:155,45:155)=NaN;
        
        % calculate bounds of image from centroid
        b = [iFr(i,2)-99 iFr(i,3)-99 iFr(i,2)+100 iFr(i,3)+100];
        
        % if frame is padded with zeros at the edge, crop the image
        inbounds = ~any(b(1:2)<1) & b(3)<camRes(2) & b(4)<camRes(1);       
        if ~inbounds

            if b(1)<1
                tmp_frame(:,(b(1):b(3)<1))=[];
                b(1)=1;
            end
            if b(2)<1
                tmp_frame((b(2):b(4)<1),:)=[];
                b(2)=1;
            end
            if b(3)>camRes(2)
                tmp_frame(:,(b(1):b(3)>camRes(2)))=[];
                b(3)=camRes(2);
            end
            if b(4)>camRes(1)
                tmp_frame((b(2):b(4)>camRes(1)),:)=[];
                b(4)=camRes(1);
            end
            
        end

        
        % extract matching region of the reference image
        refsub = ref(b(2):b(4),b(1):b(3));
        replace = isnan(refsub) & tmp_frame~=0;
        refsub(replace) = tmp_frame(replace);
        ref(b(2):b(4),b(1):b(3)) = refsub;
        end
        
    end
    
    delete(hwb);
    
    bg_lum = nanmedian(ref(:));
    segmentationOptions.bg_lum = bg_lum;
    ref = uint8(ref);
    ref(ref==0)=bg_lum;
    ref2 = uint8(ones(camRes+refpad*2).*median(double(ref(:))));
    ref2(refpad+1:refpad+camRes(1),refpad+1:refpad+camRes(2)) = ref;
    ref = ref2;
    save([fDir '\' fLabel '_ref.mat'],'ref');
    clearvars ref2 refsub tmp_frame b iFr
    
end

    
%ref = imgaussfilt(ref,parameters.refSigma);
    
    
%%  Area normalization and (possibly) bodyThreshold finding

    maxArea = 8000;
    
    % Grab random frames from input movie
    idx = randi([startImage,nFrames],[parameters.areaNormalizationNumber,1]);
    %basisImage = imresize(basisImage,[150 150]);
    basisSize = sum(basisImage(:)>0);
    s = size(basisImage);
    currentImageSet = uint8(zeros(s(1),s(2),parameters.areaNormalizationNumber));
    
    parfor i=1:length(idx)
        currentImageSet(:,:,i) = read(vidObj,idx(i));
    end
    
    if mean(currentImageSet(:)) > 100
        currentImageSet = imcomplement(currentImageSet);
    end
    %%
    % Automatically find body threshold if bodyThreshold is set to -1
    if bodyThreshold < 0   
        T = zeros(length(idx),parameters.areaNormalizationNumber);
        parfor i=1:length(idx)
            b = [cenDat(idx(i),1)-99 cenDat(idx(i),2)-99 cenDat(idx(i),1)+100 cenDat(idx(i),2)+100] + refpad;
            refsub = ref(b(2):b(4),b(1):b(3));
            [testImage,mask] = segmentDiffIm(currentImageSet(:,:,i),refsub,...
                5,.05,imageThreshold,[],[],minArea,true);
            if sum(mask(:)) > minArea && sum(mask(:)) < maxArea
                II = testImage(testImage>0);
                T(i) = autoFindThreshold_gmm(II,3);
            end
        end
        
        T = T(T>0);
        bodyThreshold = quantile(T,.5); 
        disp(['bodyThresh = ' num2str(bodyThreshold)]);
        parameters.bodyThreshold = bodyThreshold;
         
    end
    
    
    if parameters.asymThreshold < 0
        parameters.asymThreshold = parameters.bodyThreshold;
        asymThreshold = parameters.asymThreshold;
    end
            
    % Calculate the median number of pixel above the threshold
    imageSizes = zeros(size(idx));
    for j = 1:parameters.areaNormalizationNumber
        a = currentImageSet(:,:,j);
        imageSizes(j) = sum(a(:)>bodyThreshold);
    end
    imageSize = median(imageSizes);
    areanorm = sqrt(basisSize/imageSize);
    

    
 %%       
    if isempty(image_path)
        image_path   = input('Image Path = ?:  ', 's');
    end
    
    [status,~]=unix(['ls ' image_path]);
    if status == 1
        unix(['mkdir ' image_path]);
    end
    
    if ~segmentationOff
        referenceImage = segmentImage_combo(basisImage,dilateSize,...
            cannyParameter,imageThreshold,[],[],minArea,true);
    else
        referenceImage = basisImage;
    end
    
    % Define range for flipping detector
    [ii,~] = find(referenceImage > 0);
    minRangeValue = min(ii) - rangeExtension;
    maxRangeValue = max(ii) + rangeExtension;
    
    segmentationOptions.referenceImage = referenceImage;
    segmentationOptions.minRangeValue = minRangeValue;
    segmentationOptions.maxRangeValue = maxRangeValue;
    
    %define groupings
    imageVals = startImage:finalImage;
    numImages = length(imageVals);
    minNumPer = floor(numImages / numProcessors+1e-20);
    remainder = mod(numImages,numProcessors);
    count = 1;
    
    % Assign frame numbers to be aligned by each processor
    groupings = cell(numProcessors,1);
    for i=1:numProcessors
        if i <= remainder
            groupings{i} = imageVals(count:(count+minNumPer));
            count = count + minNumPer + 1;
        else
            groupings{i} = imageVals(count:(count+minNumPer-1));
            count = count + minNumPer;
        end
    end

    
    
    % Write Out Grouping Start and Finish indices
    groupidx = zeros(length(groupings),2);
    for i = 1:length(groupings)
        groupidx(i,1) = groupings{i}(1);
        groupidx(i,2) = groupings{i}(end);
    end
    
    
    %initialize new avi files
    alignmentFiles = cell(numProcessors,1);
    fDigits = ceil(log10(numProcessors+1e-10));
    for i=1:numProcessors
        qq = num2str(i);
        qq = [repmat('0',1,fDigits - length(qq)) qq];
        alignmentFiles{i} = VideoWriter([image_path '/' fLabel 'file' qq '.avi'],'Grayscale AVI');
    end
    
    x1s = zeros(numProcessors,1);
    y1s = zeros(numProcessors,1);
    angle1s = zeros(numProcessors,1);
    area1s = zeros(numProcessors,1);
    svdskip1s = zeros(numProcessors,1);
    
    currentPhis = zeros(numProcessors,1);
    
    
    %initialize First Images
    
    images = cell(numProcessors,1);

    for j=1:numProcessors
        
        i = groupings{j}(1);
        
        % If image is RGB, extract red channel only
        originalImage = read(vidObj,i);
        if length(size(originalImage)) == 3
            originalImage = originalImage(:,:,1);
        end
        % Segment the image if segmentation is ON
        if ~segmentationOff
            b = [cenDat(i,1)-99 cenDat(i,2)-99 cenDat(i,1)+100 cenDat(i,2)+100] + refpad;
            refsub = ref(b(2):b(4),b(1):b(3));
            [imageOut,mask2,bg] = segmentDiffIm(originalImage,refsub,dilateSize,cannyParameter,...
                imageThreshold,[],[],minArea,true);
            imageOut = rescaleImage(imageOut,areanorm);
        else
            imageOut = originalImage;
            imageOut = rescaleImage(imageOut,areanorm);
        end
        
        imageOut2 = imageOut;
        imageOut2(imageOut2 < asymThreshold) = 0;
        
        % If any non-zero pixels remain
        if max(imageOut2(:)) ~= 0
            
            % rotate the image to the initial Phi
            [angle1s(j),x1s(j),y1s(j),~,~,image] = ...
                alignTwoImages(referenceImage,imageOut2,initialPhi,spacing,pixelTol,false,imageOut);
            
            % Set pixels below minRange or above maxRange to zero
            s = size(image);
            b = image;
            b(b > asymThreshold  + bg) = 0;
            if minRangeValue > 1
                b(1:minRangeValue-1,:) = 0;
            end
            if maxRangeValue < length(referenceImage(:,1))
                b(maxRangeValue+1:end,:) = 0;
            end
            
            % Create fractional histogram of pixels along x-dimension
            q = sum(b) ./ sum(b(:));
            % Subtract x centroid from symmetry dividing line
            asymValue = symLine - sum(q.*(1:s(1)));
                        
            % If asymValue is negative, correct rotational degeneracy
            if asymValue < 0
                
                initialPhi = mod(initialPhi+180,360);
                
                [tempAngle,tempX,tempY,~,~,tempImage] = ...
                    alignTwoImages(referenceImage,imageOut2,initialPhi,spacing,pixelTol,false,imageOut);
                
                b = tempImage;
                b(b > asymThreshold  + bg) = 0;
                if minRangeValue > 1
                    b(1:minRangeValue-1,:) = 0;
                end
                if maxRangeValue < length(referenceImage(:,1))
                    b(maxRangeValue+1:end,:) = 0;
                end
                
                q = sum(b>0)./sum(b(:)>0);
                asymValue2 = symLine - sum(q.*(1:s(1)));
                
                if asymValue2 > asymValue
                    angle1s(j) = tempAngle;
                    x1s(j) = tempX;
                    y1s(j) = tempY;
                    image = tempImage;
                end
                
            end
            
            area1s(j) = sum(mask2(:));
            currentPhis(j) = angle1s(j);
            images{j} = image;
            svdskip1s(j) = 0;
            
        else
            
            area1s(j) = sum(mask2(:));
            currentPhis(j) = initialPhi;
            angle1s(j) = initialPhi;
            svdskip1s(j) = 1;
            image = uint8(zeros(size(imageOut2)));
            images{j} = image;
            
        end
        
    end
    
     %% Load or create new basis image

 %{
if exist([fDir '\' fLabel '_basis.png'],'file')
    
        newBasis = imread([fDir '\' fLabel '_basis.png']);
        
else
    
    tic
    Xs_temp = cell(1,1);
    Ys_temp = cell(1,1);
    Angles_temp = cell(1,1);
    Areas_temp = cell(1,1);
    svdskips_temp = cell(1,1);
    idx = randi([startImage,nFrames],[500,1]);
    i=1;
    
    [newBasis] = ...
    getBasisImage(idx,currentPhis(i),...
    segmentationOptions,nDigits,file_path,alignmentFiles{i},readout,i,...
    asymThreshold,area1s(i),vidObj,[],areanorm,images{i},ref,cenDat,refpad);
    imwrite(newBasis,[fDir '\' fLabel '_basis.png']);
    
end
 %}

segmentationOptions.referenceImage = imread([segmentationOptions.mapdir ...
    '\image_processing\basisImage.tiff']);   
%% Align the images
    
    fprintf(1,'Aligning Images\n');
    
    tic
    Xs_temp = cell(numProcessors,1);
    Ys_temp = cell(numProcessors,1);
    Angles_temp = cell(numProcessors,1);
    Areas_temp = cell(numProcessors,1);
    svdskips_temp = cell(numProcessors,1);
%%    
    parfor i=1:numProcessors
%%        
        [Xs_temp{i},Ys_temp{i},Angles_temp{i},Areas_temp{i},svdskips_temp{i}] = ...
            align_subroutine_parallel_avi_ZW(groupings{i},currentPhis(i),...
            segmentationOptions,nDigits,file_path,alignmentFiles{i},readout,i,...
            asymThreshold,area1s(i),vidObj,[],areanorm,images{i},ref,cenDat,refpad);
%%        
        
        
        Xs_temp{i}(1) = x1s(i);
        Ys_temp{i}(1) = y1s(i);
        Areas_temp{i}(1) = area1s(i);
        Angles_temp{i}(1) = angle1s(i);
        svdskips_temp{i}(1) = svdskip1s(i);
        
        close(alignmentFiles{i});   
        
    end
    
    
    Xs = combineCells(Xs_temp);
    Ys = combineCells(Ys_temp);
    angles = combineCells(Angles_temp);
    areas = combineCells(Areas_temp);
    svdskips = combineCells(svdskips_temp);
    
        
    
    x = abs(diff(unwrap(angles.*pi/180).*180/pi));
    framesToCheck = find(x > 90) + 1;
    svdskipped = find(svdskips == 1);
    

    


