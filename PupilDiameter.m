%% Measure pupil diameter for JPEG
%This script will find the pupil diameter and position in a video of the
%mouse eye. Modified from a version written by Matt McGinley.

%Basic outline to script
%1) Load jpeg of eye in for-loop
%2) Canny edge detection on original image.
%3) Fit ellipse to edge pixels
%4) Find area, centroid and Major and Minor axis of pupil ellipse
%5) Store variables in MySQL using DataJoint


%% Set path and load jpeg 

% Read jpeg images (Garrett's case)
% path for Garrett's data: /home/jantine/newnas/Garrett/Pupil Videos/010316/Mouse882/1/010316_Mouse...CompressedROIs/010316_Mouse...Eyes/010316_EyeFrames/2/... 

yourFolder = '/home/jantine/newnas/Garrett/Pupil Videos/010316/Mouse882/1/010316_Mouse882_1_CompressedROIs/010316_Mouse882_1_Eyes/010316_EyeFrames/set';
addpath(yourFolder);



% convert jpeg files into concatenated .mat file
% jpegFiles = dir('*.jpeg');
% numfiles = length(jpegFiles);
% mydata = cell(1, numfiles);
% 
% for k = 1:numfiles
%     mydata{k} = imread(jpegFiles(k).name);
% end
% 
% save('video2.mat', 'mydata');


% convert jpeg files into concatenated .mat file
% img_str = {sprintf('010316_Mouse882_1_Eye_%04d.jpeg\n', 20:29999)};
% img_names = regexp(img_str{:}(1:end-1), '\n', 'split');
% filename = 'video.mat';
% save(filename, 'img_names');
% load('video.mat')

%% Select Region of Interest (ROI) - maybe not really necessary in Garrett's case
I = imread('010316_Mouse882_1_Eye_29.jpeg');
%I = imread('010316_Mouse882_1_Eye_28757.jpeg');
I = rgb2gray(I);

% BW change to black and white only
level = 0.1; %level in range [0,1]
BW = im2bw(I, level);
figure;
imshow(BW)

%% Edge Detection (Canny)
BW1 = edge(BW, 'Canny'); 
figure;
imshow(BW1); title('Edge in BW1');
title('Canny Filter');

% Enlarge figure to full screen
%set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);


% get rid of small blobs
BW1 = bwareaopen(BW1, 220);

% get the boundary and outline over original image
boundaries = bwboundaries(BW1);
x = boundaries{1}(:, 2);
y = boundaries{1}(:, 1);
% display original image
imshow(BW1); title('Outline over Original Image');
% plot boundaries over image
hold on;
plot(x, y, 'g-', 'LineWidth', 2);
hold off

maxCoor = [min(x), max(y)];
minCoor = [max(x), min(y)];

% D1 = bwdist(BW1, 'euclidean');
% RGB1 = repmat(mat2gray(D1), [1 1 3]);
% figure
% imshow(RGB1), title('Euclidean')
% 
% RGB1 = rgb2gray(RGB1);
% BW2 = edge(RGB1, 'Canny');
% figure;
% imshow(BW1)
% title('Canny Filter, selection');

%[BW1, thresh] = edge(I, 'Canny'); 
%thresh % first element specifies lower threshold, below which all edges are disregarded; second element the higher threshold above which al edge pixels are preserved


%% Find Center and Area
% Centroid is the center of mass of the region
% Area that returns a scalar that specifies the actual number of pixels
% Eccentricity of the elipse indicate whether it is more a circle (0) or a
%   line segment (1) - validate whether Major/MinorAxisLength are correctly
%   used as the area will be transformed to a circle
% MajorAxisLength is a scalar with the length in pixels of the major axis
%   of the ellipse that has the same normalized second central moments
% MinorAxisLength is the same but then for the minor axis
% Orientation is a scalar that gives the angle between the x-axis and the
%   major axis of the ellipse
% see help regionprops for more information

% BW2_temp = regionprops(BW1, 'Image');
% BW2 = [BW2_temp.Image];

s = regionprops(BW1, 'Centroid', 'Area', 'Eccentricity', ...
    'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Image');
centroids = cat(1, s.Centroid);
majAxis = cat(1, s.MajorAxisLength);
minAxis = cat(1, s.MinorAxisLength);
area = cat(1, s.Area);
orientation = cat(1,s.Orientation);


% plot centroid and boundaries in original image
imshow(I); title('Outline Original with RegionProps');
hold on
plot(x, y, 'g-', 'LineWidth', 1);
plot(centroids(:,1), centroids(:,2), 'r*');
  
% hlen = minAxis/2;
% xCenter = centroids(1);
% yCenter = centroids(2);
% cosOrient = cosd(orientation);
% sinOrient = sind(orientation);
% xcoords = xCenter + hlen .* [cosOrient -cosOrient];
% ycoords = yCenter + hlen .* [-sinOrient sinOrient];
% line(xcoords, ycoords);

line(maxCoor, minCoor)

hold off



% get diameter and plot in original image - TO SHOW THAT DATA IS
% TRANSFORMED

diameters = mean([s.MajorAxisLength s.MinorAxisLength], 2);
radii = diameters/2;
% plot the circles
viscircles(centroids,radii);



% get area that returns a scalar that specifies the actual 

% get centers and radii of the circles
% stats = regionprops('table', BW1, 'Centroid', ...
%     'MajorAxisLength', 'MinorAxisLength');
% centers = stats.Centroid;
% diameters = mean([stats.MajorAxisLength stats.MinorAxisLength], 2);
% radii = diameters/2;
% diameter = diameters;
% % plot the circles
% viscircles(centers,radii);
% hold off



%% Get diameter measurement
% We can use the MajorAxisLength and MinorAxisLength to compare across the
% pupils to look for changes. As getting the diameter/radii circles
% transforms your data much. Maybe we need to get more data points of the
% pupil ellipse to compare. 

% and compare area
 
%% Adjusted Dov script for jpeg (not working properly

% %adjusted Dov script
%         [row1, col1] = find(I);
%         pupilCenter_estimate = [median(row1),median(col1)]; %center
%         pupilRadius_estimate = sqrt(length(row1)/pi);   %radius (from A=pi*r^2)
%                 
%         radiusFudge = 1.4;  %fudge factor for estimated radius; default 1.4
%         edgeFudge = 1;   %default 1 (lower means more edges will be found)
%         threshPupilDist = 3;  %edge pixels must be at least this close to dark pixels; default 3
%         minEdgeSize = 15; %edge pixels must be continguous with at least this many other pixels; default 15
%         %find edges in image with canny image detection
%         [~, threshold] = edge(I, 'canny');
%         im_edges = edge(I,'canny', threshold * edgeFudge);
%         %delete edge pixels which are too far from center
%         [row2, col2] = find(im_edges);
%         for i=1:length(row2)
%             %distance from estimated center
%             dXY = pupilCenter_estimate - [row2(i), col2(i)];
%             d = sqrt(sum(dXY.*dXY));
%             %delete if too far
%             if d>pupilRadius_estimate*radiusFudge
%                 im_edges(row2(i),col2(i)) = false;
%             else
%                 %delete edge pixels which are too far from dark pixels
%                 [~,d] = dsearchn([row1,col1],[row2(i),col2(i)]);
%                 if d>threshPupilDist
%                     im_edges(row2(i),col2(i)) = false;
%                 end
%             end
%         end
%         
%         %delete edge pixels which are not part of large group
%         im_edges2 = bwlabel(im_edges); %label ROIs by unique number
%         %number of ROIs in edges matrix
%         nRois = max(im_edges2(:));
%         %clear edge matrix and add back in if ROI is big
%         im_edges = false(size(im_edges));
%         for iROI = 1:nRois
%             roiSize = sum(sum(im_edges2==iROI));
%             if roiSize >= minEdgeSize
%                 im_edges(im_edges2==iROI) = true;
%             end
%         end
%         
%         %only try ellipse fitting if there are sufficient edges pixels
%         if sum(im_edges(:))>minEdgeSize
%             %fit ellipse
%             [row3,col3]=find(im_edges);
%             pupilEllipse = fit_ellipse(row3,col3);
%             %if ellipse found, use short axis as pupil diameter
%             if ~isempty(pupilEllipse.long_axis) && pupilEllipse.long_axis~=0
%                 %record pupil diameter
%                 pupilD(I) = pupilEllipse.long_axis;
%                 %record pupil position
%                 pupilXY(1,I) = pupilEllipse.Y0_in;
%                 pupilXY(2,I) = pupilEllipse.X0_in;
%             end
%         end
%         
%         %plot
%         imshow(I);
%         %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%         %pupil outline
%         if ~isempty(pupilEllipse)
%             ellipse(pupilEllipse.b,pupilEllipse.a,pupilEllipse.phi,pupilEllipse.Y0_in,pupilEllipse.X0_in,'r');
%             %pupil center
%             hold on
%             scatter(pupilEllipse.Y0_in,pupilEllipse.X0_in,'r');
%             hold off
%         end

