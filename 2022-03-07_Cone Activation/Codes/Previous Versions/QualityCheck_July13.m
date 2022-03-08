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
unstabilizedVideo=VideoReader([PathName,filesep,videoName(1:end-11),'.avi']);
fNom=unstabilizedVideo.FrameRate*unstabilizedVideo.Duration;
FindPRL_July13
crossLoc_raw=nanmean(alllocs_OneStrip_raw);
load([videoName(1:end-11),'_meanrem_960_hz_34.mat'])
alllocs_OneStrip_all=alllocs_OneStrip_raw;
crossFlagOne=find(~isnan(alllocs_OneStrip_raw(:,1))==1);
crossFlagStartOne=crossFlagOne(1);
crossFlagEndOne=crossFlagOne(end);
alllocs_OneStrip_all(crossFlagOne(1):crossFlagOne(end),:)=frameshifts_strips((crossFlagOne(1)-1)*33+17:33:(crossFlagOne(end)-1)*33+17,:)+crossLoc_raw+[100,100];
FindPRL_February13
if ~isempty(setdiff(crossFlagStartOne:crossFlagEndOne,crossFlagOne))
    cprintf('red',['Cross_finder couldn''t find cross in frame ',num2str(setdiff(crossFlagStartOne:crossFlagEndOne,crossFlagOne)),'\n'])
%     beep
end

figure;
subplot(1,3,1)
plot(crossFlagStartOne:crossFlagEndOne,-ceil(frameshifts_strips((crossFlagStartOne-1)*33+16:33:(crossFlagEndOne-1)*33+16,1))+crossLoc_raw(1),'o')
hold on
plot(crossFlagStartOne:crossFlagEndOne,alllocs_OneStrip(crossFlagStartOne:crossFlagEndOne,1)-100,'.')

subplot(1,3,2)
plot(crossFlagStartOne:crossFlagEndOne,-round(frameshifts_strips((crossFlagStartOne-1)*33+16:33:(crossFlagEndOne-1)*33+16,1))+crossLoc_raw(1),'o')
hold on
plot(crossFlagStartOne:crossFlagEndOne,alllocs_OneStrip(crossFlagStartOne:crossFlagEndOne,1)-100,'.')

subplot(1,3,3)
plot(crossFlagStartOne:crossFlagEndOne,-floor(frameshifts_strips((crossFlagStartOne-1)*33+16:33:(crossFlagEndOne-1)*33+16,1))+crossLoc_raw(1),'o')
hold on
plot(crossFlagStartOne:crossFlagEndOne,alllocs_OneStrip(crossFlagStartOne:crossFlagEndOne,1)-100,'.')


-ceil(frameshifts_strips(16:33:end,:))+alllocs_OneStrip_raw