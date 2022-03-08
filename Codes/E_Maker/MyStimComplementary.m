close all
clear
clc
load('BAK1008L_2020_07_03_11_23_38_AOSLO.mat')
tGuess = -0.6;
tGuessSd = 0.3;
beta=3; delta=0.01; gamma=0.5;
pThreshold = (1-gamma)/3*2+gamma;
q.normalizePdf=1;
correct=1;

trial=1;
ntrials=20;
files=dir(fullfile(cd, '*.avi'));
files={files.name};
times_all_orientations=floor(ntrials/8);
remaining_trials=ntrials-(floor(ntrials/8)*8);
% E_orientations=[];
dotSize=CFG.dotsize;
dotSeparation=CFG.dotseparation;

% for trial_num=1:times_all_orientations
%     E_orientations=[E_orientations; ((ceil(randperm(8)/2))*90-90)'];
% end
% E_orientations=[E_orientations; ((ceil(randperm(8,remaining_trials)/2))*90-90)'];
% E_orientations=rawData(:,3);
% q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma);
% intensity_log=QuestQuantile(q);	% Recommended by Pelli (1987), and still our favorite.
% intensity=10^intensity_log*10;    %=intensity/gap size or offset in pixel
offset=-rawData(:,2);
% offsety=0;
% offset=CFG.first;
VLnew=MycreateStimulus(offset(1),dotSize,dotSeparation,0,0);
myVideo=VideoReader(cell2mat(files(1)));
CurrentFrame=read(myVideo,10);
sz(1,1:2)=size(VLnew);
if mod(sz(1,1),2)==0
    VLnew=padarray(VLnew,[1 1],1,'post');
end
figure;
subplot(1,2,1);imshow(VLnew,[])
subplot(1,2,2);imshow(CurrentFrame,[])
sz(1,3:4)=size(VLnew);
allVLC(:,:,1)=1-VLnew;
for trial=2:ntrials
    clear VLnew
    correct=rawData(trial-1,5);
    
    VLnew=MycreateStimulus(offset(trial),dotSize,dotSeparation,0,0);
    myVideo=VideoReader(cell2mat(files(2*trial-1)));
    CurrentFrame=read(myVideo,10);
    sz(trial,1:2)=size(VLnew);
    if mod(sz(trial,2),2)==0
        VLnew=padarray(VLnew,[0 1],1,'post');
    end
    sz(trial,3:4)=size(VLnew);
    figure;
    subplot(1,2,1);imshow(VLnew,[])
    subplot(1,2,2);imshow(CurrentFrame,[])
    VLpad=padarray((1-VLnew),[0 (size(allVLC(:,:,1),2)-size(VLnew,2))/2],0,'both');
    allVLC(:,:,trial)=VLpad;
end
save('VLsC.mat','allVLC');
