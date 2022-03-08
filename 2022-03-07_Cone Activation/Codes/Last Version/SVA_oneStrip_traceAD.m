% This code makes cone activation map for the input AOSLO video and saves
% the results in a structure
% Asieh Daneshi
% last modification: 05.11.2021
tic
close all
clearvars -except alllocs_OneStrip alllocs_OneStrip_unstab alllocs_OneStrip_30frames newStimAll crossFlagOne crossFlagStartOne crossFlagEndOne videoName videoNames PathName imName CDC20_loc boxposition conelocs sumNorm_Jenny alllocs_OneStrip_backup matName1 picMat_areas picMat_activities ConeLoc_actLogs alllocs_OneVideo alllocs_OneVideo_unstab alllocs_OneVideo_30frames       
% clc
warning('off','all')
%%
% keptFrames=setdiff(1:30,Outliers);   % removing outliers. 30 is the number of frames
keptFrames=find(~isnan(alllocs_OneStrip(:,1)));
clc
% fprintf('Please select the video you aim to analyze. (offline stabilized video)\n');
% [videoName,PathName]=uigetfile('*stabilized_del.avi','MultiSelect', 'off');
% [videoName,PathName]=uigetfile('*_stab_A_stabilized.avi','MultiSelect', 'off');
myVideo=VideoReader([PathName,filesep,videoName]);
fNom=myVideo.FrameRate*myVideo.Duration;    % Number of frames in the video
% sumNormMaker
matchCDC
toc
% PRLX=WCentroid20(1);   % PRL location
% PRLY=WCentroid20(2);
% clc
%% ========================================================================
Im=sumNorm_Jenny;
ref=im2double(read(myVideo,crossFlagStartOne));
[~,varargout]=corr2d(Im,ref);
[rref,cref]=find(ref);
% clc
% Removing conelocs outside boxposition -----------------------------------
PRLX=356;   % PRL location
PRLY=356;  % PRL location
conelocs(conelocs(:,1)<boxposition(1),:)=[];
conelocs(conelocs(:,1)>boxposition(1)+boxposition(3),:)=[];
conelocs(conelocs(:,2)<boxposition(2),:)=[];
conelocs(conelocs(:,2)>boxposition(2)+boxposition(4),:)=[];
% Shifting the cone locations to the upper right corner of the image (for
% later use on single frames, not sumnorm) --------------------------------
conelocsN=conelocs;
% conelocsN(:,1)=conelocsN(:,1)-boxposition(1)-100;
% conelocsN(:,2)=conelocsN(:,2)-boxposition(2)-100;
%% ========================================================================
% =========================================================================
% Initiating variables ----------------------------------------------------
ZF=1;
cropSize=710*ZF+1;
FullField=512;  % In pixel of imaging grid
HalfFF=floor(FullField/2);    % Middle point of FullField
% ap_field=FullField;	% The field size, in 'pixels' (typically 512)
ap_field=cropSize;	% The field size, in 'pixels' (typically 512)
% zernike_pupil is the pupil size for the wavefront testing (usually==psf_pupil):
[psf_pupil,zernike_pupil]=deal(7);	% The pupil size for testing
% "deal" simply matches up the input and output lists. For example,
% [a,b,c,...]=deal(5), mean a=5,b=5,c=5
Zoom=10;	% draw PSF much smoother by "zooming in"
% ZF=200;   % additional Zoom Factor (used for cone delivery)
PixPerDegree=600;   % It means each arcmin is 10 pixels
% field_size is the field size, in ARCMIN; (60 arcmin per degree, so 72 for a 1.2deg field):
field_size=(FullField/PixPerDegree*60);    % In 'arcmin'
% field_size=(FullField/PixPerDegree*60)/Zoom;    % In 'arcmin'
% lambda=.788;    % The stimulus wavelength;
lambda=.840;    % The stimulus wavelength;
diff_limited=1; % if diff_limited==1, set all Zernike coefficients to zero; else, input a given .zer from HSWS
defocus=0.03;  % in diopter, 0.01?
% -------------------------------------------------------------------------
% Building PSF
myPSF=GeneratePSF(ap_field,psf_pupil,zernike_pupil,field_size,lambda,diff_limited,defocus);


%% ========================================================================
% ZF=1;  % zoom factor
[v,c]=voronoin(conelocs(:,1:2)*ZF);    % N-D Voronoi diagram

A=zeros(length(c),1) ;  % Matrix that will contain area of each voronoi cell

% computing the area of each voronoi cell ---------------------------------
for i=1:length(c)
    v1=v(c{i},1) ;
    v2=v(c{i},2) ;
    A(i)=polyarea(v1,v2);
end
%--------------------------------------------------------------------------
%% ========================================================================
ref=imresize(im2double(read(myVideo,1)),ZF);
sz=size(ref);
S=zeros(sz);

% cropSize is the size of that part of image, around PRL that we are
% focusing at. Please only assign odd numbers to cropSize -----------------
cropSize=710*ZF+1;
conelocIm=zeros(cropSize,cropSize);

% croppedFrames=zeros(cropSize,cropSize);



CurrentFrame=imresize(im2double(read(myVideo,1)),ZF);
conelocsNN=conelocs;
% ---------------------------------------------------------------------
% croppedFrames=CurrentFrame(PRLX*ZF-(cropSize-1)/2:PRLX*ZF+(cropSize-1)/2,PRLY*ZF-(cropSize-1)/2:PRLY*ZF+(cropSize-1)/2);
clear croppedCL
% croppedCL=conelocsNN((conelocsNN(:,1)>=PRLX*ZF-(cropSize-1)/2) & (conelocsNN(:,1)<=PRLX*ZF+(cropSize-1)/2) & (conelocsNN(:,2)>=PRLY*ZF-(cropSize-1)/2) & (conelocsNN(:,2)<=PRLY*ZF+(cropSize-1)/2),:);
croppedCL=conelocsNN;
croppedConeLocs=croppedCL;
% croppedConeLocs(:,1,1)=croppedConeLocs(:,1,1)-PRLX*ZF+(cropSize-1)/2;
% croppedConeLocs(:,2,1)=croppedConeLocs(:,2,1)-PRLY*ZF+(cropSize-1)/2;
[croppedVoronoi_v,croppedVoronoi_c]=voronoin(croppedConeLocs(:,1:2));    % N-D Voronoi diagram
croppedVoronoi_c(end-1)=croppedVoronoi_c(end-2);
croppedVoronoi_c(end)=croppedVoronoi_c(end-1);
for n=1:length(croppedVoronoi_c)
    croppedVoronoi_c{n,1}(length(croppedVoronoi_c{n,1})+1)=croppedVoronoi_c{n,1}(1);
end
currentA=zeros(length(croppedVoronoi_c),1);	% Matrix that will contain area of each voronoi cell
for i=1:length(croppedVoronoi_c)
    v1=croppedVoronoi_v(croppedVoronoi_c{i},1);
    v2=croppedVoronoi_v(croppedVoronoi_c{i},2);
    currentA(i)=polyarea(v1,v2);
    if ~isnan(currentA(i)) && (currentA(i)<=500*(ZF)^2)        
        %&&(currentA(i)<=350*(ZF)^2)
        currentD(i)=maxDiameter(croppedVoronoi_c{i},croppedVoronoi_v);   % Cone diameter
        SpotSize=floor(2*currentD(i));  % Ask Niklas!
        Aperture=currentD(i)*0.48;   % proportion of inner segement diameter (ISD) that functions a light collection aperture (FWHM) from MacLeod et al, 1992; set range as desired
        cA=Aperture./2.25482;
        HalfSS=floor(SpotSize/2);   % Half of SpotSize
        spot=single(fspecial('gaussian',SpotSize,cA));	%generate Gaussian (SpotSize=ConeSize pixels wide)
        spot=spot./max(spot(:));	%normalize Gaussian
        sz_conelocIm=size(conelocIm);
        % In this loop we go for each cone location in the selected
        % part of the current frame -----------------------------------
        croppedConeLocs1=ceil(croppedConeLocs(i,1,1));
        croppedConeLocs2=ceil(croppedConeLocs(i,2,1));
        if(croppedConeLocs1>0 && croppedConeLocs2>0)    % remove meaningless cone locations
            tempConeIm=zeros(cropSize); % Make a whole black image
            % Whiten the only one pixel in the position of the current cone
            tempConeIm(ceil(croppedConeLocs(i,2,1)),ceil(croppedConeLocs(i,1,1)))=1;
            % Convolve the image of the current cone with a PSF
            filtSpot=conv2(tempConeIm,spot,'same');
            C=zeros(size(filtSpot));
            [rt,ct]=find(filtSpot~=0);
            % Find out which parts of the convolved PSF is inside the current voronoi cell
            inVoronoi=inpolygon(ct,rt,croppedVoronoi_v(croppedVoronoi_c{i},1),croppedVoronoi_v(croppedVoronoi_c{i},2));
            filtSpotEs=zeros(sz_conelocIm(1:2));
            rIn=rt(inVoronoi);
            cIn=ct(inVoronoi);
            for j2=1:length(rIn)
                filtSpotEs(rIn(j2),cIn(j2))=filtSpot(rIn(j2),cIn(j2));
            end
            conelocIm=conelocIm+filtSpotEs;
        end
    end
end

allA(1:length(currentA),1)=currentA;

% preparing an avi file to write videos ===================================
vWrite=VideoWriter([videoName(1:end-4),'_stim_one.avi']);
vWrite.FrameRate=30;
open(vWrite);
% load([videoName(1:end-4),'_stabilized.mat'])      % loading mat file containing all the stimuli
myStim=1-newStimAll;

ConeLoc_actLog = nan(length(croppedConeLocs),fNom);
croppedConeLocs_rnd = round(croppedConeLocs);
croppedConeLocs_rnd(croppedConeLocs_rnd == 0) = 1;

[r_out,c_out] = find(alllocs_OneStrip(:,1)<min(conelocs(:,1))+size(myStim,1)/2 | alllocs_OneStrip(:,1)>max(conelocs(:,1))-size(myStim,1)/2 | alllocs_OneStrip(:,2)<min(conelocs(:,2))+size(myStim,2)/2 | alllocs_OneStrip(:,2)>max(conelocs(:,2))-size(myStim,2)/2);
alllocs_OneStrip_rnd = round(alllocs_OneStrip);

a2=0;
a3=0;
for a1=1:fNom
        CurrentFrame=imresize(im2double(read(myVideo,a1)),ZF);
        if ~isempty(find(CurrentFrame))&& ismember(a1,crossFlagOne) && ~ismember(a1,r_out)
%         if ~isempty(find(CurrentFrame))&& ismember(a1,crossFlagOne) && alllocs_OneStrip(a1,2)*ZF-size(myStim,1)/2*ZF>boxposition(1) && alllocs_OneStrip(a1,2)*ZF+size(myStim,1)/2*ZF<boxposition(3) && alllocs_OneStrip(a1,1)*ZF-size(myStim,2)/2*ZF>boxposition(2) && alllocs_OneStrip(a1,1)*ZF+size(myStim,2)/2*ZF<boxposition(4)
            a3=a3+1;
            stimIm_prep=zeros(cropSize);
            stimIm_prep(alllocs_OneStrip_rnd(a1,2)*ZF-floor(size(myStim,1)/2)*ZF:alllocs_OneStrip_rnd(a1,2)*ZF+floor(size(myStim,1)/2)*ZF,alllocs_OneStrip_rnd(a1,1)*ZF-floor(size(myStim,2)/2)*ZF:alllocs_OneStrip_rnd(a1,1)*ZF+floor(size(myStim,2)/2)*ZF) = myStim(:,:,a3);  % For square stimulus
%             stimIm_prep(alllocs_OneStrip(a1,2)*ZF-size(myStim,1)/2*ZF:alllocs_OneStrip(a1,2)*ZF+size(myStim,1)/2*ZF,alllocs_OneStrip(a1,1)*ZF-size(myStim,2)/2*ZF:alllocs_OneStrip(a1,1)*ZF+size(myStim,2)/2*ZF)=imresize(myStim(:,:,a3),[size(myStim,1)*ZF+1,size(myStim,2)*ZF+1]);  % For square stimulus
%             stimIm = stimIm_prep;
            stimIm = imrotate(stimIm_prep(:,:,1),Rotation, 'crop');
            %             stimIm_translated=stimIm(PRLY*ZF-(cropSize-1)/2:PRLY*ZF+(cropSize-1)/2,PRLX*ZF-(cropSize-1)/2:PRLX*ZF+(cropSize-1)/2);
            if ~isempty(find(stimIm~=0))
                a2=a2+1;
                indStim(a2)=a1;
                % =========================================================================
                % Convolving the stimulus image with the PSF ------------------------------
                filtStim=CONVOLVE(stimIm,myPSF);
                filtStim(filtStim<0.01)=0;
                OutputIm=conelocIm.*(filtStim);
                [rt,ct]=find(filtStim~=0);
                croppedConeLocstemp=croppedConeLocs((croppedConeLocs(:,1)>=min(ct(:))) & (croppedConeLocs(:,1)<=max(ct(:))) & (croppedConeLocs(:,2)>=min(rt(:))) & (croppedConeLocs(:,2)<=max(rt(:))),:);
                
                rr=zeros(1,size(croppedConeLocstemp,1));
                for m=1:size(croppedConeLocstemp,1)
                    rr(m)=min(find(croppedConeLocstemp(m,1)==croppedConeLocs(:,1)&croppedConeLocstemp(m,2)==croppedConeLocs(:,2)));
                end
                %% ========================================================================
                OutputImF=zeros(cropSize,cropSize);
                %         ii=1;
                for i=rr
                    v1=croppedVoronoi_v(croppedVoronoi_c{i},1);
                    v2=croppedVoronoi_v(croppedVoronoi_c{i},2);
                    currentA(i)=polyarea(v1,v2);
                    if ~isnan(currentA(i))  && (currentA(i)<=500*(ZF)^2)
%                         &&(currentA(i)<=350*(ZF)^2)
                        % In this loop we go for each cone location in the selected
                        % part of the current frame -----------------------------------
                        croppedConeLocs1=ceil(croppedConeLocs(i,1));
                        croppedConeLocs2=ceil(croppedConeLocs(i,2));
                        if(croppedConeLocs1>0 && croppedConeLocs2>0)    % remove meaningless cone locations
                            tempConeIm=zeros(cropSize); % Make a whole black image
                            % Whiten the only one pixel in the position of the current cone
                            tempConeIm(ceil(croppedConeLocs(i,1,1)),ceil(croppedConeLocs(i,2,1)))=1;
                            
                            
                            % Find out which parts of the convolved PSF is inside the current voronoi cell
                            inVoronoiF=inpolygon(ct,rt,croppedVoronoi_v(croppedVoronoi_c{i},1),croppedVoronoi_v(croppedVoronoi_c{i},2));
                            filtSpotF=zeros(sz_conelocIm(1:2));
                            
                            rIn=rt(inVoronoiF);
                            cIn=ct(inVoronoiF);
                            for j2=1:length(rIn)
                                filtSpotF(rIn(j2),cIn(j2))=OutputIm(rIn(j2),cIn(j2));
                            end
                            [rz,cz]=find(filtSpotF~=0);
                            mySum=sum(filtSpotF(:));
                            filtSpotEs_sum=zeros(sz_conelocIm(1:2));
                            for j3=1:length(rIn)
                                filtSpotEs_sum(rIn(j3),cIn(j3))=mySum;
                            end
                            OutputImF=OutputImF+filtSpotEs_sum;
                        end
                    end
                end
                %% ========================================================================
                OutImF=OutputImF;
                OutputImF=OutputImF/max(OutImF(:));
                picMat_area(1:cropSize,1:cropSize)=OutputImF;
                picMat_activity(1:cropSize,1:cropSize)=OutputIm;
                
                picMat_areas.(['V' videoName(end-13:end-11)]).(['F' num2str(a1,'%02d')]) = picMat_area;
                picMat_activities.(['V' videoName(end-13:end-11)]).(['F' num2str(a1,'%02d')]) = picMat_activity;
                
                for cLoop = 1:length(croppedConeLocs)                    
                    ConeLoc_actLog(cLoop,a1) = picMat_area(croppedConeLocs_rnd(cLoop,1),croppedConeLocs_rnd(cLoop,2));
                end
                
%                 save(['picMat_areaOne_',videoName(end-13:end-11),'_',num2str(a1),'.mat'],'picMat_area');
%                 save(['picMat_activityOne_',videoName(end-13:end-11),'_',num2str(a1),'.mat'],'picMat_activity');
                fig1=figure;imshow(1-OutputImF,[])
                myFrame=getframe;
                writeVideo(vWrite,myFrame)
                close(fig1)
            end
    end
end


% ConeLoc_actLog(ConeLoc_actMap < 0.001) = NaN;
ConeLoc_actLogs.(['V' videoName(end-13:end-11)]) = ConeLoc_actLog;
save('ConeLoc_actLogs.mat','ConeLoc_actLogs');

Voronoi_areas =allA;
save('Voronoi_areas.mat','Voronoi_areas');

save('picMat_areas.mat','picMat_areas');
save('picMat_activities.mat','picMat_activities');

alllocs_OneVideo.(['V' videoName(end-13:end-11)]) = alllocs_OneStrip;
save('alllocs_OneVideo.mat','alllocs_OneVideo');

alllocs_OneVideo_unstab.(['V' videoName(end-13:end-11)]) = alllocs_OneStrip_unstab;
save('alllocs_OneVideo_unstab.mat','alllocs_OneVideo_unstab');

alllocs_OneVideo_30frames.(['V' videoName(end-13:end-11)]) = alllocs_OneStrip_30frames;
save('alllocs_OneVideo_30frames.mat','alllocs_OneVideo_30frames');

save('croppedVoronoi_c.mat','croppedVoronoi_c');
save('croppedVoronoi_v.mat','croppedVoronoi_v');
save('keptFramesOne.mat','keptFrames');
save('crossFlagOne.mat','crossFlagOne');
% save(['alllocs_OneStrip_',videoName(end-13:end-11),'.mat'],'alllocs_OneStrip');
% save('alllocs_OneStrip.mat','alllocs_OneStrip');
save('keptFramesOne.mat','keptFrames');
