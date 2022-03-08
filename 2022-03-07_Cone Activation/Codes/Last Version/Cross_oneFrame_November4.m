function [x,y,crossFlag]=Cross_oneFrame_November4(videoIn,cross_num,frameN,SearchThreshold)
%This program serves as a core function to find a stimulus cross in a
%stabilized AOSLO video. Crosses for IR, red, and green stimuli have
%different pixel intensities. To search for a specific cross, use the
%following flags (the program defaults to IR):
%
%IR:    0;
%Red:   1;
%Green: 2;
%
%Input:
%(image,cross_num,SearchThreshold)
%
%Output:
%[x,y]
%

switch nargin
    case 4
        % Do nothing
        SearchThresh=SearchThreshold;
    case 3
        SearchThresh = 0.7;
    case 2
        cross_num = 10; %default is to search for IR cross
        SearchThresh = 0.7;
end
filter = zeros(11); filter(:,6)=1/11; %vertical cross
filter(6,:) = 1/11;


% [x,y,fr] = deal(zeros(endframe,1));

if cross_num == 0          % IR, if gain ~= 0
    pixel = 255/255;
elseif cross_num == 10     % IR, if gain == 0
    pixel = 254/255;
elseif cross_num == 1      % Red (digi1), if gain ~= 0
    pixel = 253/255;
elseif cross_num == 11      % Red (digi1), if gain == 0
    pixel = 252/255;
elseif cross_num == 2      % Green (digi2), if gain ~= 0
    pixel = 251/255;
elseif cross_num == 12      % Green (digi2), if gain == 0
    pixel = 250/255;
end

crossFlag=1;
currentframe=im2double(read(videoIn,frameN));
[row,col] = find(currentframe(:,:,1)==pixel);
searchframe = zeros(712,712);
for rownum = 1:size(row,1)
    searchframe(row(rownum),col(rownum))=1;
end
meanframe = normxcorr2_mex(filter,searchframe,'same');
mx=max(meanframe(:));
if mx>SearchThresh
    [ ~ , row_max , column_max ] = max2D_RS(meanframe);
    x=column_max;
    y=row_max;
else
    %         fprintf('No Cross in the whole video! Sorry\n')
    crossFlag=0;
    x=[];
    y=[];
    % delete zeros
    % ZeroEntries = x==0;
    % x(ZeroEntries) = [];
    % y(ZeroEntries) = [];
end


