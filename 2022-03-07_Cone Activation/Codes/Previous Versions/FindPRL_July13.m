% clc
%% ========================================================================
% [videoName,PathName]=uigetfile('*stabilized.avi','MultiSelect', 'off');
% fprintf('Please select the video you aim to find the PRL in it. (It should be online stabilized)\n');
% [videoName,PathName]=uigetfile('*stab_A.avi','MultiSelect', 'off');
% fprintf('Please select the stabilized video (we need it for background image).\n');
% [videoName_stab,PathName_stab]=uigetfile('.avi','MultiSelect', 'off');
% myVideo=VideoReader([PathName,filesep,videoName]);
% fNom=myVideo.FrameRate*myVideo.Duration;    % Number of frames in the video
% =========================================================================
% first we find the cross locations on the frames containing crosses ======
shiftVar=zeros(fNom,3); % Contains the shift values for each frame compared to Ref
alllocs_OneStrip_raw=NaN(fNom,2);
for a1=1:fNom
    shiftVar(a1,1)=a1;
    [x,y,crossFlag]=Cross_oneFrame_November4(unstabilizedVideo,10,a1);
    if crossFlag==1
        shiftVar(a1,2:3)=[x,y];
        alllocs_OneStrip_raw(a1,2)=shiftVar(a1,3);
        alllocs_OneStrip_raw(a1,1)=shiftVar(a1,2);
    end
end
%% ========================================================================
% finding cross locations in the remaining frames -------------------------
x=nanmedian(alllocs_OneStrip_raw(:,1));
y=nanmedian(alllocs_OneStrip_raw(:,2));
% crossFlagOne=find(~isnan(alllocs_OneStrip_raw(:,1))==1);
% crossFlagStartOne=crossFlagOne(1);
% crossFlagEndOne=crossFlagOne(end);