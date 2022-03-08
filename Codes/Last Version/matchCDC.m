% close all
clc

% fprintf('Please select the CDC location file.\n');
% [matName,PathName]=uigetfile('BAK*150.mat','MultiSelect', 'off');
% load([PathName,filesep,matName])
% fprintf('Please select the cone location file.\n');
% [matName,PathName]=uigetfile('BAK*JLR.mat','MultiSelect', 'off');
% load([PathName,filesep,matName])
% fprintf('Please select the sumnorm image(Jenny''s).\n');
% [imName,PathName]=uigetfile('BAK*.tiff','MultiSelect', 'off');
% sumNorm_Jenny=imread([PathName,filesep,imName]);

% fprintf('Please select your sumnorm image.\n');
% [imName,PathName]=uigetfile('BAK*stabilized.tif','MultiSelect', 'off');
% sumNorm_Es=imread([PathName,filesep,imName]);
sumNorm_Es=imread([videoName(1:end-11),'.tiff']);
figure;
imshow(sumNorm_Jenny,[])
hold on
plot(conelocs(:,1),conelocs(:,2),'y.')
plot(WCentroid20(1),WCentroid20(2),'g*')

sumNorm_Jenny_crop=sumNorm_Jenny(WCentroid20(1)-200:WCentroid20(1)+200,WCentroid20(2)-200:WCentroid20(2)+200);

sumNorm_Es_crop=sumNorm_Es(356-200:356+200,356-200:356+200);
Nxcorr_Es=normxcorr2(sumNorm_Es_crop,sumNorm_Jenny);
[xEs,yEs]=find(Nxcorr_Es==max(Nxcorr_Es(:)));        % cross correlation
CDC_Es=[yEs,xEs]-[401,401]+[201,201];
figure;
imshow(sumNorm_Es,[])
hold on
plot(conelocs(:,1)-WCentroid20(1)+CDC_Es(1),conelocs(:,2)-WCentroid20(2)+CDC_Es(2),'y.')
% plot(conelocs(:,1),conelocs(:,2),'y.')
conelocs(:,1)=conelocs(:,1)-WCentroid20(1)+CDC_Es(1);
conelocs(:,2)=conelocs(:,2)-WCentroid20(2)+CDC_Es(2);
plot(CDC_Es(1),CDC_Es(2),'g+')
for n=1:length(alllocs_OneStrip)
    if ~isnan(alllocs_OneStrip(n,1))
        alllocs_OneStrip(n,1)=alllocs_OneStrip(n,1)-WCentroid20(1)+CDC_Es(1);
        alllocs_OneStrip(n,2)=alllocs_OneStrip(n,2)-WCentroid20(2)+CDC_Es(2);
        plot(alllocs_OneStrip(n,1),alllocs_OneStrip(n,2),'ws')
    end
end