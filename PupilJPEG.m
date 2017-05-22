%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% README
%
% author: Jantine Broek
% e-mail: jantine.broek@yale.edu
% date: May 2017
% for: McCormick lab, Yale University, New Haven, USA
%
%
% Pre-req:
%            fit_ellipse  -  from https://www.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse
%                            Conic Ellipse representation = a*x^2+b*x*y+c*y^2+d*x+e*y+f=0
%                            (Tilt/orientation for the ellipse occurs when the term x*y exists (i.e. b ~= 0)) 
%                            EDIT: made changes in this function to create
%                            plots.
%
%            natsortfiles -  from http://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort
%
%            ellipse      -  from https://www.mathworks.com/matlabcentral/fileexchange/289-ellipse-m.
%                            Only used to draw ellipse in figure.
%
% Goal:
%           (1) Using edge detection, the pupil is found and the coordinates of the
%           pupil location will be used to estimate the function for (2) ellipse
%           fitting. With the formula for ellipse fitting, the (3) regionprops will
%           be used to compare the dynamics of the pupil fluctuation.
%
%
% Input:    
%           .JPEG - these are files of individual pupil images obtained with Virtual Dub 
%                   software (http://www.virtualdub.org/download.html)             
%
%
% Output:   
%           longAxis   - size of the long axis of the ellipse
%           shortAxis  - size of the short axis of the ellipse
%           pupilXY    - pupil location
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
cd /home/jantine/newnas/Garrett/'Pupil Videos'/010316/Mouse882/1/010316_Mouse882_1_CompressedROIs/010316_Mouse882_1_Eyes/010316_EyeFrames/2
yourFolder = '/home/jantine/newnas/Garrett/Pupil Videos/010316/Mouse882/1/010316_Mouse882_1_CompressedROIs/010316_Mouse882_1_Eyes/010316_EyeFrames/2';
% cd /Volumes/'Seagate Backup Plus Drive'/Pupil_Garrett/set
% yourFolder = '/Volumes/Seagate Backup Plus Drive/Pupil_Garrett/set';
addpath(yourFolder);

% place files in natural order 
contents = dir('*.jpeg'); 
n = natsortfiles({contents.name});
[~,ndx] = natsortfiles({contents.name}); 
contents = contents(ndx);

% load data in directory of MatLab path
filename = cell(1,numel(n));
for j = 1:numel(n)
  filename(j) = {contents(j).name};
  I = imread(filename{j});
  
  % display loop progression
  if mod(j,1000)==0
      disp([num2str(j) ' / ' num2str(size(n)) ' image ']);
  end
  
end

%% Select eye ROI

% average images
picEyeAvg = mean(I,3);

% choose region of interest around eye
figure
    uiwait(msgbox('Select region of interest. Double click to close loop. '));
    maskEye = roipoly(picEyeAvg./255);
close all

%% Choose threshold for pupil value 

% threshold
threshDark = 25;

% check threshold value
figure
uiwait(msgbox('Choose threshDark value so that most pixels of pupil are showing.'));
for iTest = 1:15
    randFrame = ceil(rand*size(filename,2));
    disp(['random frame #' num2str(randFrame)])
    
    randImg = imread(filename{randFrame});
    randImg = randImg(:,:,1);
    
    randImg_threshed = randImg<threshDark;
    imshow(randImg_threshed.*maskEye)
    pause
end
close all

  
%% Edge detection and Ellipse fitting

% eye may be closed if very few dark pixels found (default = 200)
threshNdarkPix = 200; 

% fudge factor for estimated radius (default = 1.4)
radiusFudge = 1.4;  

% canny edge detection threshold fudge factor (default 1). Lower means more edges will be found
edgeFudge = 1;  

% edge pixels must be at least this close to dark pixels (default 3)
threshPupilDist = 3;  

% edge pixels must be continguous with at least this many other pixels (default 15)
minEdgeSize = 10;     


tic
for k = 1:numel(n)
  filename(k) = {contents(k).name};
  I = imread(filename{k});
  
  % create title for images
  [pathstr, name, ext] = fileparts(filename{k});
  ix=strfind(name,'_');
  t = name(ix(4)+1:end);
  
  
  % display loop progression
  if mod(k,500)==0
      disp([num2str(k) ' / ' num2str(size(n)) ' image ']);
  end

  
  % define image frame and threshold
    im_threshed = (I<threshDark) .* maskEye;
    
  % estimate pupil center and size based of number of pixels
    [row1, col1] = find(im_threshed);
    
    
  % check to see if eye is open
    if length(row1) > threshNdarkPix
        
        % pupil center
        pupilCenter_estimate = [median(row1), median(col1)]; 
        
        % pupil radius (from A=pi*r^2)
        pupilRadius_estimate = sqrt(length(row1)/pi);   
        
        
        % delete pixels which are too far away from pupil
        for l = 1:length(row1)
            
            % distance from estimated center
            dXY = pupilCenter_estimate - [row1(l), col1(l)];
            d = sqrt(sum(dXY .* dXY));
            
            % delete pixels if too far
            if d > pupilRadius_estimate * radiusFudge
                im_threshed(row1(l), col1(l)) = false;
            end
            
        end
        
        
   % canny edge detection
     im = rgb2gray(I);
     [~, threshold] = edge(im, 'canny');
     im_edges = edge(im, 'canny', threshold * edgeFudge) .* maskEye;
     
   % delete edge pixels which are too far from center
     [row2, col2] = find(im_edges);
     for m = 1:length(row2)
         
        % distance from estimated center
          dXY = pupilCenter_estimate - [row2(m), col2(m)];
          d = sqrt(sum(dXY .* dXY));
          
        % delete edge pixels which are too far from dark pixels
          [~,d] = dsearchn([row1, col1], [row2(m), col2(m)]);
          if d > threshPupilDist
             im_edges(row2(m),col2(m)) = false;
          end

     end       
        
    
    % delete edge pixels which are not part of large group
    % label ROIs by unique number
      im_edges2 = bwlabel(im_edges); 
      
    % number of ROIs in edges matrix
      nRois = max(im_edges2(:));
      
    % clear edge matrix and add back in if ROI is big
      im_edges = false(size(im_edges));
      for iROI = 1:nRois
          roiSize = sum(sum(im_edges2 == iROI));
          
          if roiSize >= minEdgeSize
             im_edges(im_edges2==iROI) = true;
          end
          
      end
      
    % ellipse fitting if there are sufficient edges pixels
        if sum(im_edges(:)) > minEdgeSize
            
           % fit ellipse
             [row3, col3] = find(im_edges);
             pupilEllipse = fit_ellipse(row3, col3);
             
            % if ellipse found, use long axis as pupil diameter
            if ~isempty(pupilEllipse.long_axis) && pupilEllipse.long_axis~=0
                
                % record pupil axis
                shortAxis(k) = pupilEllipse.short_axis;
                longAxis(k) = pupilEllipse.long_axis; 
                
                % record pupil position
                pupilXY(1,k) = pupilEllipse.Y0_in;
                pupilXY(2,k) = pupilEllipse.X0_in;
            else
                % no ellipse found but parabola or hyperbola
                error_fit(k) = filename(k); 
            end
            
        end
        
    end         

        
        
%   % plot pic with fitted ellipse
%     figure;
%     imshow(im); title(sprintf('Pupil Ellipse 2, image %s', t));
%     
%     % resize figure
%     set(gcf, 'Position', [150, 150, 400, 400]); 
%     
%     if ~isempty(pupilEllipse)
%             
%         ellipse(pupilEllipse.b,pupilEllipse.a,pupilEllipse.phi,pupilEllipse.Y0_in,pupilEllipse.X0_in,'r');
%             
%         % pupil center
%         hold on 
%         scatter(pupilEllipse.Y0_in,pupilEllipse.X0_in,'r'); 
%             
%         % pupil axis
%         cos_phi = pupilEllipse.cos_phi;
%         sin_phi = pupilEllipse.sin_phi;
% 
%         % rotation matrix to rotate the axes with respect to an angle phi
%         R = [ cos_phi sin_phi; -sin_phi cos_phi ];
% 
%         % the axes
%         ver_line        = [ [pupilEllipse.X0 pupilEllipse.X0]; pupilEllipse.Y0+pupilEllipse.b*[-1 1] ];
%         horz_line       = [ pupilEllipse.X0+pupilEllipse.a*[-1 1]; [pupilEllipse.Y0 pupilEllipse.Y0] ];
%         new_ver_line    = R * ver_line;
%         new_horz_line   = R * horz_line;
% 
%         % draw major axis
%         plot( new_ver_line(2,:),new_ver_line(1,:),'b' ); 
%         % draw minor axis
%         plot( new_horz_line(2,:),new_horz_line(1,:),'g' ); 
% 
%         hold off
%      end

end
toc

%% Extract data tuples for DataJoint
% extract session date, session number, filename for x-axis, movie number (first number of filenr), frame
% number (rest of numbers in filenr), animal id
token = strtok(filename,'.');
D = regexp(token, '_', 'split');
D = vertcat(D{:});

% filenr
filenr_temp = D(:,5);
G = sprintf('%s*', filenr_temp{:}); % change cell to double
filenr = sscanf(G, '%f*');
filenr = filenr';

% session date
session_date_temp = D(:,1);
session_date = datetime(session_date_temp, 'InputFormat', 'ddMMyy', 'Format', 'yyyy-MM-dd');
session_date = string(session_date); % as string
session_date = char(session_date); % as char temporarily for DataJoint
session_date = session_date';

% session number
session_number_temp = D(:,3);
F = sprintf('%s*', session_number_temp{:}); % change cell to double
session_number = sscanf(F, '%f*');


% movie number and frame number
for m = 1:numel(filenr)
    G = num2str(filenr(m));
    movie_number(m) = str2double(G(1));
    frame_number(m) = str2double(G(2:end));
end

% extract animal_id column
animal_id_temp = D(:,2);
animal_id_temp = animal_id_temp(:,1);
animal_id_temp = sprintf('%s*', animal_id_temp{:});
animal_id = regexp(animal_id_temp, '\d*', 'Match');
animal_id = str2double(animal_id');

%% Plot axis
figure;
hold on

subplot(2,1,1); 
plot(shortAxis, '-.ob')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
legend('shortAxis');

subplot(2,1,2); 
plot(longAxis, '-.or')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
legend('longAxis');
hold off

%% Collect data
data = [filenr; shortAxis; longAxis; pupilXY]; 
data(data == 0) = NaN;
data = data';

% save data
fname = sprintf('pupilAxis_Mouse%d.mat', animal_id(1:1));
save(fname, 'data');

% find NaN
[row, col] = find(isnan(data));
nrNaN = unique(row);
dataNaN = data(nrNaN,:); % all is NaN
% check NaN images
filename_NaN = filename(:,nrNaN)';
figure; 
for ii = 1:numel(filename_NaN)
   INaN = imread(filename_NaN{ii});
   imshow(filename_NaN{ii})
   pause;
end


%% Outlier removal

% compute the mean value
shortAxis_mean = nanmean(shortAxis);
longAxis_mean = nanmean(longAxis);
pupilXY_mean = nanmean(pupilXY,2);

% compute the absolute difference
shortAxis_absdiff = abs(shortAxis - shortAxis_mean);
longAxis_absdiff = abs(longAxis - longAxis_mean);
pupilXY_absdiff = abs(pupilXY - pupilXY_mean);

% compute the median of the absolute difference
shortAxis_mad = median(shortAxis_absdiff);
longAxis_mad = median(longAxis_absdiff);
pupilXY_mad = median(pupilXY_absdiff);

% outliers if the absolute difference is moe than some factor times the mad
% value
% shortAxis
sensitivityFactor = 4; %change this value for more critical detection
thresholdValue_sA = sensitivityFactor * shortAxis_mad;
outlierIndexes_sA = abs(shortAxis_absdiff) > thresholdValue_sA;
% longAxis
thresholdValue_lA = sensitivityFactor * longAxis_mad;
outlierIndexes_lA = abs(longAxis_absdiff) > thresholdValue_lA;
% pupilXY
thresholdValue_XY = sensitivityFactor * pupilXY_mad;
outlierIndexes_XY = abs(pupilXY_absdiff) > thresholdValue_XY;

% extract outlier values
shortAxis_outliers = shortAxis(outlierIndexes_sA);
longAxis_outliers = longAxis(outlierIndexes_lA);
pupilXY_outliers = pupilXY(outlierIndexes_XY);

% % extract non-outlier values
% shortAxis_nonoutliers = shortAxis(~outlierIndexes_sA);
% longAxis_nonoutliers = longAxis(~outlierIndexes_lA);
% pupilXY_nonoutliers = pupilXY(~outlierIndexes_XY);

% remove outliers from original dataset and replace with NaN
shortAxis(outlierIndexes_sA) = nan;
longAxis(outlierIndexes_lA) = nan;
%pupilXY(outlierIndexes_XY) = nan;

%% Interpolate NaN values
shortAxis_inter = inpaintn(shortAxis);
longAxis_inter = inpaintn(longAxis);

% plot axis with interpolated values
figure;
plot(data(:,3), '-or')
hold on
plot(longAxis_inter, '-og')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
title('Original and interpolated pupil long axis')
legend('Original', 'LP filtered');
hold off


%% TO DO
% check pupilXY outliers






%% Interpolate missing values
data_inter = data;
data_inter = inpaintn(data_inter);
% check that NaN are interpolated
% [row1, ~] = find(isnan(data_inter));
data_inter_ex = data_inter(nrNaN,:);
shortAxis_NaN = data_inter(:,2);
longAxis_NaN = data_inter(:,3);


% data_inter2 = repnan(data_inter2);
% % check that NaN are interpolated
% % [row1, ~] = find(isnan(data_inter));
% data_inter2_ex = data_inter2(nrNaN,:);


% calculate mean pupil position and delta position over time
pupilXY_inter(1,:) = inpaintn(pupilXY(1,:));
pupilXY_inter(2,:) = inpaintn(pupilXY(2,:));
pupilXY_mean = nanmean(pupilXY_inter,2);

% difference of interpolated NaNs with mean (delta)
pupilXYDel = nan(1,length(pupilXY_inter));
for ij = 1:length(pupilXY_inter)
    if ~isnan(pupilXY_inter(1,ij)) 
        ydel = pupilXY_mean(1) - pupilXY_inter(1,ij);
        xdel = pupilXY_mean(2) - pupilXY_inter(2,ij);
        pupilXYDel(ij) = sqrt(ydel^2 + xdel^2);
    end
end


% locate average pupil center location
figure;
imshow(mean(I,3),[0 255])
hold on
scatter(pupilXY_mean(1),pupilXY_mean(2));
hold off

figure;
plot(pupilXYDel+40,'Linewidth',2);
title('Pupil Diameter and Eye Position');
xlabel('Time (1/10 sec)');
ylabel('Pupil Diameter (pixels)');
legend('Eye position');

% plot axis with interpolated values
figure;
plot(longAxis, '-or')
hold on
plot(longAxis_NaN2, '-og')
xticks(1:1:numel(filenr));
set(gca,'XTickLabel',filenr)
xtickangle(45)
xlabel('filenr');
ylabel('length axis (px)');
title('Original and interpolated pupil long axis')
legend('Original', 'LP filtered');
hold off



% calculate mean pupil position and delta position over time
pupilXY_mean = nanmean(pupilXY_inter,2);

%delta position over time
pupilXYDel = nan(1,length(pupilXY_inter));
for ij=1:length(pupilXY_inter)
    if ~isnan(pupilXY_inter(1,ij)) %avoid using values we deleted manually
        ydel = pupilXY_mean(1) - pupilXY_inter(1,ij);
        xdel = pupilXY_mean(2) - pupilXY_inter(2,ij);
        pupilXYDel(ij) = sqrt(ydel^2 + xdel^2);
    end
end

% imshow(mean(video,3),[0 255])
% hold on
% scatter(pupilXY_mean(1),pupilXY_mean(2));
% hold off
% 
% figure
% plot(pupilD_inter,'Linewidth',2);
% hold on
% plot(pupilXYDel+40,'Linewidth',2);
% hold off
% title('Pupil Diameter and Eye Position');
% xlabel('Time (1/10 sec)');
% ylabel('Pupil Diamter (pixels)');
% legend('Pupil diameter','Eye position');

% Low-pass filter pupil trace
%frequency cuttoff for low-pass filter
freqCutoff = 4;
sampleInterval = mode(diff(videoTS));
[b,a]=besself(4,2*pi*freqCutoff);
[bd,ad] = bilinear(b,a,1/sampleInterval);
pupilD_smooth = filtfilt(bd,ad,pupilD_inter-mean(pupilD_inter));
pupilD_smooth = pupilD_smooth + mean(pupilD_inter);

figure
plot(pupilD,'Linewidth',2)
hold on
plot(pupilD_smooth,'r','Linewidth',2)
title('Original and filtered pupil diameter measurements')
legend('Original', 'LP filtered');
hold off

correlogram = xcorr(pupilD_smooth,pupilD_inter,50,'coeff');
%plot(-50:1:50,correlogram)


%% Add to DataJoint
% addpath /media/jantine/Data/04_DataJoint/2PE/schemas
% 
% % collect data for schema preprocess.EyeROI
% tuple.animal_id = animal_id;
% tuple.session_date = session_date; % should be in string format "2017-05-12"
% tuple.session_number = session_number;
% tuple.movie_number = movie_number';
% tuple.frame_number = frame_number';
% tuple.edge_area = edgeArea';
% tuple.short_axis = shortAxis';
% tuple.long_axis = longAxis';
% 
% % order structure as DataJoint keys
% a = preprocess.EyeROI;
% fields = a.header.names;
% tuple = orderfields(tuple, fields);
% 
% insert(preprocess.EyeROI, dj.struct.fromFields(tuple))
