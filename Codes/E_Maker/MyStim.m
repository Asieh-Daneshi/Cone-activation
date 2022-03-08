close all
clear
clc
rawData=[1,2,0,0,1,1;2,2.50000000000000,90,90,1,2;3,3,180,180,1,3;4,3.50000000000000,270,270,1,4;5,4,0,0,1,5;6,5,90,90,1,6;7,6,180,180,1,7;8,7,270,270,1,8];
save('rawData.mat','rawData');
tGuess=-0.25;
tGuessSd=0.2;
beta=3; delta=0.01; gamma=0.25;
pThreshold=(1-gamma)/3*2+gamma;
q.normalizePdf=1;
correct=1;

trial=1;
ntrials=8;
files=dir(fullfile(cd, '*.avi'));
files={files.name};
times_all_orientations=floor(ntrials/8);
remaining_trials=ntrials-(floor(ntrials/8)*8);
E_orientations=[];

% for trial_num=1:times_all_orientations
%     E_orientations=[E_orientations; ((ceil(randperm(8)/2))*90-90)'];
% end
% E_orientations=[E_orientations; ((ceil(randperm(8,remaining_trials)/2))*90-90)'];
E_orientations=rawData(:,3);
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma);
intensity_log=QuestQuantile(q);	% Recommended by Pelli (1987), and still our favorite.
intensity=10^intensity_log*10;    %=intensity/gap size or offset in pixel
offset=rawData(:,2);
offsety=0;
% offset=CFG.first;
Enew=MycreateStimulus(offset(1),0,0,E_orientations(1),0);
% myVideo=VideoReader(cell2mat(files(1)));
% CurrentFrame=read(myVideo,10);
sz(1,1:2)=size(Enew);
if mod(sz(1,1),2)==0
    Enew=padarray(Enew,[1 1],1,'post');
end
% figure;
% subplot(1,2,1);imshow(Enew,[])
% subplot(1,2,2);imshow(CurrentFrame,[])
sz(1,3:4)=size(Enew);

allE(:,:,1)=padarray((1-Enew),[(71-length(Enew))/2 (71-length(Enew))/2],0,'both');
for trial=2:ntrials
    clear Enew
    correct=rawData(trial-1,5);
%     [offset,alt]=newAlt(intensity, correct,alt, trial, E_orientations);
    
    Enew=MycreateStimulus(offset(trial),0,0,E_orientations(trial),0);
%     myVideo=VideoReader(cell2mat(files(trial)));
%     CurrentFrame=read(myVideo,10);
    sz(trial,1:2)=size(Enew);
    if mod(sz(trial,1),2)==0
        Enew=padarray(Enew,[1 1],'replicate','post');
    end
    sz(trial,3:4)=size(Enew);
    figure;
    imshow(Enew,[])
%     subplot(1,2,2);imshow(CurrentFrame,[])
    Epad=padarray((1-Enew),[(71-length(Enew))/2 (71-length(Enew))/2],0,'both');
    allE(:,:,trial)=Epad;
end
save('Es.mat','allE');
