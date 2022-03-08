
clear
clc
fprintf('Please select the video you aim to analyze. (online stabilized video)\n');
[videoName,PathName]=uigetfile('*stabilized.avi','MultiSelect', 'off');
myVideo=VideoReader([PathName,filesep,videoName]);
fNom=myVideo.FrameRate*myVideo.Duration;  
for a=1:fNom
    I(:,:,a)=im2double(read(myVideo,a));
end
[sx,sy,sz]=size(I);
for a1=1:sx
    for a2=1:sy
        sumNorm(a1,a2)=sum(I(a1,a2,:))/length(nonzeros(I(a1,a2,:)));
    end
end
figure;imshow(sumNorm,[])
imwrite(sumNorm,[videoName(1:end-4),'.tiff'])