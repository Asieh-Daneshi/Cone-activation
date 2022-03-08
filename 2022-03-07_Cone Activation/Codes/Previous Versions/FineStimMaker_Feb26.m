close all
clear a1 n ref shiftVar x y
%% ========================================================================

load('Es.mat')
% myStim=padarray(1-allE(:,:,str2double(videoName(end-13:end-11))),[0 5],'replicate','post');
myStim=1-allE(:,:,str2double(videoName(end-13:end-11)));
count=0;
isnotnan=find(~isnan(alllocs_OneStrip(:,1))==1);
for b=1:length(find(~isnan(alllocs_OneStrip(:,1))==1))
    a=isnotnan(b);
    count=count+1;
    Im=im2double(read(myVideo,a));
    ImCrop=Im(alllocs_OneStrip(a,2)-(size(myStim,2)-1)/2:alllocs_OneStrip(a,2)+(size(myStim,2)-1)/2,:);
    ImCrop(ImCrop<=0.004)=0;
    mask=ImCrop;
    mask(ImCrop~=0)=1;
    for b=1:size(myStim,2)
        [~,c]=find(mask(b,:));
        [(max(c)+min(c))/2,alllocs_OneStrip(a,1)];
        
        nr=1;
        nc=size(myStim,2);
        Nr=ifftshift((-fix(nr/2):ceil(nr/2)-1));
        Nc=ifftshift((-fix(nc/2):ceil(nc/2)-1));
        [Nc,Nr]=meshgrid(Nc,Nr);
        Greg=fft2(myStim(b,:)).*exp(1i*2*pi*((0*Nr/nr)+((alllocs_OneStrip(a,1)-(max(c)+min(c))/2)*Nc/nc)));
%         Greg=Greg.*exp(-1i*0);
        
        newStim(b,:)=ifft2(Greg);
        newStim(newStim<=0.1)=0;
    end
%     newStimAll(:,:,count)=newStim;
    newStimAll(:,:,count)=myStim;
%     figure;
%     subplot(3,1,1);imshow(ImCrop,[])
%     subplot(3,1,2);imshow(mask,[])
%     subplot(3,1,3);imshow(newStim,[])
end
save([videoName(1:end-4),'_stabilized.mat'],'newStimAll')
%%
delete *stabilized.avi
vWrite=VideoWriter([videoName(1:end-4),'_stabilized.avi'],'Grayscale AVI');    % start a video to record jitter-free frames
vWrite.FrameRate=30;
open(vWrite);


vRead=VideoReader([PathName,filesep,videoName]);
fNom=vRead.FrameRate*vRead.Duration;
for b=1:fNom
    currentFrame=read(vRead,b);
    writeVideo(vWrite,uint8(currentFrame))
end
close(vWrite)