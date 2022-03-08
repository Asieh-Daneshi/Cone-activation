close all
clear
clc
%%
current_folder=pwd;
files=dir(fullfile(current_folder,'BAK1041R*.csv'));
NameList={files.name};
NameList=sort(NameList);
for a1=1:length(NameList)
    currentName=cell2mat(NameList(a1));
    CSVfile=csvread(currentName);
    Xs=CSVfile(1:30,4:2:69);
    Ys=CSVfile(1:30,5:2:69);
    allShifts=NaN(size(Xs,1)*size(Xs,2),5);
    for a2=1:size(Xs,1)
        for a3=1:size(Xs,2)
            allShifts((a2-1)*33+1:a2*33,1)=a2;  % frame number
            allShifts((a2-1)*33+a3,2)=a3;   % strip number
            if mod(a3,33)==1 || mod(a3,33)==0
                timing((a2-1)*33+a3,1)=1000/(512*30)*8;
            else
                timing((a2-1)*33+a3,1)=1000/(512*30)*16;
            end    
        end
        allShifts(1,3)=0;
        allShifts(1,4)=timing(1,1)/2;
        allShifts(1,5)=timing(1,1);
        for a4=2:length(timing)
            allShifts(a4,3)=allShifts(a4-1,5); % timing (start of the strip)
            allShifts(a4,5)=timing(a4,1)+allShifts(a4-1,5); % timing (finish of the strip)
            allShifts(a4,4)=mean([allShifts(a4,5),allShifts(a4,3)]); % timing
        end
    allShifts(:,6)=reshape(Xs',990,1);  % x position
    allShifts(:,7)=reshape(Ys',990,1);  % y position
    end
    save([currentName(1:end-4),'_allShifts.mat'],'allShifts');
    csvwrite(['Traces_',currentName(1:end-4),'.csv'],'allShifts')
end