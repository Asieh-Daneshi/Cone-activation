% close all
clear
clc
% cone1 (miliseconds) 
tr=0.05;
td=0.3;
tp=0.4;
phi=-58.12;
% lambda_max=840
% J0=lambda/lambda_max;
J0=1;
t=0:1:1000;
J=J0*((t/1000)/tr).^3/(1+((t/1000)/tr).^3).*exp(-((t/1000)/td).^2).*cos(2*pi*(t/1000)/tp+phi);
Jp=J;
figure;plot(t,J,'r')

%% ========================================================================
load('Voronoi_areas.mat')
load('allConeActivity.mat')
AllConeActivities=AllConeActivity{1,1};    % cone activity in frames
%% ========================================================================
totalTime=1000; % mseconds, for 30*512*512 pixels
[r,c]=find(AllConeActivities~=0);

for n=1:length(r)
    a=r(n);
    ConeActivity=AllConeActivities(a,:);
    ConeActivity_ms=reshape([ConeActivity;zeros(32,30)],1,[]);

    ConeActivityNZ=find(ConeActivity_ms~=0);
    conv_ConeActivity_ms=nan(1,30*round(1000/30)+length(J)-1);
    for b=1:length(ConeActivityNZ)
        ConeActivity_ms_temp=zeros(1,30*round(1000/30)); 
        ConeActivity_ms_temp(ConeActivityNZ(b))=ConeActivity_ms(ConeActivityNZ(b));
        conv_ConeActivity_ms_new=conv(ConeActivity_ms_temp,J);
        conv_ConeActivity_ms=nanmax(conv_ConeActivity_ms,conv_ConeActivity_ms_new);
    end
    figure;
    subplot(1,3,1);bar(ConeActivity_ms)
    subplot(1,3,2);bar(conv_ConeActivity_ms)
    subplot(1,3,3);bar(conv(ConeActivity_ms,J))
end