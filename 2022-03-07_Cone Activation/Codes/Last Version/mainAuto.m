% This code is the main code for "cone activation"
% Asieh Daneshi September 2021
close all
clear
clc
%%
fprintf('Please select the CDC location file.\n');
[matName1,PathName]=uigetfile('BAK*150.mat','MultiSelect', 'off');
load([PathName,filesep,matName1])
fprintf('Please select the cone location file.\n');
[matName,PathName]=uigetfile('BAK*JLR.mat','MultiSelect', 'off');
load([PathName,filesep,matName])
fprintf('Please select the sumnorm image(Jenny''s).\n');
[imName,PathName]=uigetfile('BAK*.tiff','MultiSelect', 'off');
sumNorm_Jenny=imread([PathName,filesep,imName]);

% clearvars -except WCentroid20 boxposition conelocs sumNorm_Jenny

%% ========================================================================
fprintf('Please select all the videos you aim to analyze. (offline stabilized videos)\n');
[videoNames,PathName]=uigetfile('*_stab_A.avi','MultiSelect', 'on');
if ~iscell(videoNames)
    videoName=videoNames;
    videoName_unstab=[videoName(1:end-11),'.avi'];
    QualityCheck_July14
    FineStimMaker_July16
    SVA_oneStrip_traceAD
else
    for Nom=1:length(videoNames)
        clear alllocs_OneStrip cone_locs
        load([PathName,filesep,matName1])
        videoName=cell2mat(videoNames(Nom));
        videoName_unstab=[videoName(1:end-11),'.avi'];
        QualityCheck_July14
        FineStimMaker_July16
        SVA_oneStrip_traceAD
    end
end
