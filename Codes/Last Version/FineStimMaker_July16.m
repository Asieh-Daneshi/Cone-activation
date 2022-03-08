close all
clear a1 n ref shiftVar x y
%% ========================================================================

load('Es.mat')
% myStim=padarray(1-allE(:,:,str2double(videoName(end-13:end-11))),[0 5],'replicate','post');
myStim=1-allE(:,:,str2double(videoName(end-13:end-11)));
count=0;
isnotnan=find(~isnan(alllocs_OneStrip(:,1))==1);
for b=1:length(find(~isnan(alllocs_OneStrip(:,1))==1))
    newStimAll(:,:,b)=myStim;
end
% save([videoName(1:end-4),'_stabilized.mat'],'newStimAll')
% %%
% delete *stabilized.avi
% vWrite=VideoWriter([videoName(1:end-4),'_stabilized.avi'],'Grayscale AVI');    % start a video to record jitter-free frames
% vWrite.FrameRate=30;
% open(vWrite);
% 
% 
% vRead=VideoReader([PathName,filesep,videoName]);
% fNom=vRead.FrameRate*vRead.Duration;
% for b=1:fNom
%     currentFrame=read(vRead,b);
%     writeVideo(vWrite,uint8(currentFrame))
% end
% close(vWrite)