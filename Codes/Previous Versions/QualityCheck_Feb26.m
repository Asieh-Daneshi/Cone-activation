% close all
% clear
clc
% [videoName,PathName]=uigetfile('*stab_A.avi','MultiSelect', 'off');
myVideo=VideoReader([PathName,filesep,videoName]);
fNom=myVideo.FrameRate*myVideo.Duration;
for n=1:fNom
    ref=im2double(read(myVideo,n));
    ref(ref<=0.004)=0;
%     mask=ref;
%     mask(ref~=0)=1;
%     figure;imshow(ref,[])
%     figure;
%     subplot(1,2,1);imshow(ref,[])
%     subplot(1,2,2);imshow(mask,[])
end
%% ========================================================================
FindPRL_July14
if ~isempty(setdiff(crossFlagStartOne:crossFlagEndOne,crossFlagOne))
    cprintf('red',['Cross_finder couldn''t find cross in frame ',num2str(setdiff(crossFlagStartOne:crossFlagEndOne,crossFlagOne)),'\n'])
%     beep
end