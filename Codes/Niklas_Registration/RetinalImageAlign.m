function [ShiftX,ShiftY,Rotation,ImB_New_uint8] = RetinalImageAlign(varargin)
%% Try aligning AOSLO retinal images from different sessions
% if you are providing the input always pair wise, 2 or 3 images
% Path1, Filname1, Path2, Filname2 
% or
% Path1, Filname1, Path2, Filname2, Path3, Filname3

if isempty(varargin)
    StartPath  = 'M:\Documents\Projekte\Accommodation_RetinalStretching_dewarped\dewarped_Images(aligned_lowerRightCorner)\';
    [F1,WorkPath1] = uigetfile('*.tiff;*.tif', 'Select Reference image', StartPath);
    [F2,WorkPath2] = uigetfile('*.tiff;*.tif', 'Select Experiment image', WorkPath1);
    
else
    NumFiles = size(varargin,2)/2;
    for ax = 1:NumFiles
        eval(['WorkPath' num2str(ax) ' = varargin{ax*2-1};'])
        eval(['F' num2str(ax) ' = varargin{ax*2};'])
    end
    if NumFiles==2
        AnotherThird=1;
    else
        AnotherThird=0;
    end
end

% create Save name
if strcmp(F1(end-4:end),'.tiff')
    RefFN = F1(1:end-5);
else
    RefFN = F1(1:end-4);
end
if strcmp(F2(end-4:end),'.tiff')
    SaveFN = F2(1:end-5);
else
    SaveFN = F2(1:end-4);
end

fprintf([F2 '\non\n' F1 '\n'])

preImA = imread([WorkPath1 F1]);
preImB = imread([WorkPath2 F2]);

ImA = double(preImA(:,:,1))./255;       % fft2 etc. only work with doubles and singles
ImB = double(preImB(:,:,1))./255;

CropB = [100 125 150 200];
% CropB = [470 450 50 50];
rot_angles = -2:.2:2;
BlurSigma = [0 16 24 32];

Aligning = 1;
TryWithThirdRetinaImage = 0;

ReDo = 0;
while Aligning
    
    FullProcess = 1;
    
    % search for existing data
    EI_folder_content = dir(fullfile(WorkPath2,[SaveFN '*_AlignReady.mat']));
    if ~isempty(EI_folder_content)
        load([WorkPath2 EI_folder_content.name])
        
        % already processed?
        for ax = 1:size(AlignData,2)
            if strcmp(AlignData(ax).Reference,RefFN)
                if ~ReDo
                    FullProcess = 0;
                end
                ShiftX = AlignData(ax).ShiftX;
                ShiftY = AlignData(ax).ShiftY;
                Rotation = AlignData(ax).Rotation;
                GoodMatch = 1;
            end
        end
    end
        
        
    if FullProcess
        if TryWithThirdRetinaImage
            if AnotherThird
                [F3,WorkPath3] = uigetfile('*.tiff;*.tif', 'Select Experiment image', WorkPath2);
            end
            preImC = imread([WorkPath3 F3]);
            ImC = double(preImC(:,:,1))./255;

            [GoodMatch_1,ShiftX_1,ShiftY_1,Rotation_1] = RetinalImageAlign_core(ImA,ImC,BlurSigma,CropB,rot_angles);
            [GoodMatch_2,ShiftX_2,ShiftY_2,Rotation_2] = RetinalImageAlign_core(ImC,ImB,BlurSigma,CropB,rot_angles);

            fprintf([F3 '\non\n' F1 '\nand\n' F2 '\non\n' F3 '\n'])

            if GoodMatch_1 && GoodMatch_2
                GoodMatch = 1;

                ShiftX = ShiftX_1+ShiftX_2;
                ShiftY = ShiftY_1+ShiftY_2;
                Rotation = Rotation_1+Rotation_2;

            else
                GoodMatch = 0;
                if ~AnotherThird
                    AnotherThird = 1;
                end

            end

        else
            [GoodMatch,ShiftX,ShiftY,Rotation] = RetinalImageAlign_core(ImA,ImB,BlurSigma,CropB,rot_angles);
        end
    end



    %% apply found adjustments for alignment
    Okay = 0;
    if GoodMatch

        fprintf('Xshift: %2.1f; Yshift: %2.1f;\nRotation: %1.1f\n',...
            ShiftX,ShiftY,Rotation)

        ImB_Rot = imrotate(preImB(:,:,1),Rotation,'crop');
        deltar = -ShiftY;
        deltac = -ShiftX;
        phase = 2;
        [nr,nc]=size(ImB);
        Nr = ifftshift((-fix(nr/2):ceil(nr/2)-1));
        Nc = ifftshift((-fix(nc/2):ceil(nc/2)-1));
        [Nc,Nr] = meshgrid(Nc,Nr);
        ImB_New = ifft2(fft2(ImB_Rot).*exp(1i*2*pi*(deltar*Nr/nr+deltac*Nc/nc))).*exp(-1i*phase);

        ImB_New_double = abs(ImB_New);
        ImB_New_uint8 = uint8(round(ImB_New_double));



        %% show result
        CheckFig = figure;
        set(CheckFig, 'KeyPressFcn','uiresume');
        CheckingResults = 1;
        ImA_large = imresize(ImA(101:612,101:612),2);
        ImB_large = imresize(ImB_New_uint8(101:612,101:612),2);
        imshow(ImA_large)
        title(F1,'FontSize',12,'Color',[0 0 0],'Interpreter','none')
        xlabel('Space = next image; o = okay; n = not okay','FontSize',12,'Color',[0 0 0])

        Showing = 1;

        while CheckingResults
            uiwait(CheckFig);
            resp = get(CheckFig,'CurrentKey');
            if strcmp(resp,'o');    % okay
                CheckingResults = 0;
                Okay = 1;
                Aligning = 0;
                disp('okay')

            elseif strcmp(resp,'n');    % not okay
                CheckingResults = 0;  
                disp('not okay')
                
                [ShiftX,ShiftY,Rotation] = deal([]);
                
                if ~FullProcess
                    ReDo = 1;
                    disp('Rerunning Alignment!')
                    
                else
                    
                    if ~AnotherThird && ~TryWithThirdRetinaImage
                        TryWithThirdRetinaImage = 1;
                        AnotherThird = 0;
                        
                    elseif ~AnotherThird && TryWithThirdRetinaImage
                        AnotherThird = 1;
                        
                    else
                        % adjust parameters
                        Ans = questdlg('Can you provide another image from the same location?','Alignment failed','Yes','No','Yes');
                        
                        if strcmp(Ans,'No')
                            Aligning = 0;
                            
                        else
                            TryWithThirdRetinaImage = 1;
                        end
                    end
                end

            elseif strcmp(resp,'space');    % switch image
                if Showing == 1
                    imshow(ImB_large)
                    title(F2,'FontSize',12,'Color',[0 0 0],'Interpreter','none')
                    xlabel('Space = next image; o = okay; n = not okay','FontSize',12,'Color',[0 0 0])
                    Showing = 2;

                else
                    imshow(ImA_large)
                    title(F1,'FontSize',12,'Color',[0 0 0],'Interpreter','none')
                    xlabel('Space = next image; o = okay; n = not okay','FontSize',12,'Color',[0 0 0])
                    Showing = 1;
                end
            end
        end

        close(CheckFig)
        disp('-------------------------------------------------')
        
        
    else
        Ans = questdlg('Can you provide another image from the same location?','Alignment failed','Yes','No','Yes');
        
        if strcmp(Ans,'No')
            Aligning = 0;
            
        else
            TryWithThirdRetinaImage = 1;
        end
    end
end



%% save data
if Okay && FullProcess
    ds = datestr(now,'yyyy_mm_dd');

    imwrite(ImB_New_uint8,[WorkPath2 SaveFN '_AlignReady_' ds '.tif'])
    
    EI_folder_content = dir(fullfile(WorkPath2,[SaveFN '*_AlignReady.mat']));
    if ~isempty(EI_folder_content)
        load([WorkPath2 EI_folder_content.name])
        
        if ReDo
            for ax = 1:size(AlignData,2)
                if strcmp(AlignData(ax).Reference,RefFN)
                    AlignData(ax).Reference = RefFN;
                    AlignData(ax).ShiftX = ShiftX;
                    AlignData(ax).ShiftY = ShiftY;
                    AlignData(ax).Rotation = Rotation;
                end
            end
            
        else
            CSize = size(AlignData,2);
            AlignData(CSize+1).Reference = RefFN;
            AlignData(CSize+1).ShiftX = ShiftX;
            AlignData(CSize+1).ShiftY = ShiftY;
            AlignData(CSize+1).Rotation = Rotation;
        end
        
    else
        
        AlignData.Reference = RefFN;
        AlignData.ShiftX = ShiftX;
        AlignData.ShiftY = ShiftY;
        AlignData.Rotation = Rotation;
        
    end
    save([WorkPath2 SaveFN '_AlignReady.mat'], 'AlignData');
        
    fprintf('Data saved\n-------------------------------------------------\n')

else
    fprintf('No data was saved\n------------------------------------------\n')
end