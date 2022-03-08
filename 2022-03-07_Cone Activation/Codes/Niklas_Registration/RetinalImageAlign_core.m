function [GoodMatch,ShiftX_out,ShiftY_out,Rotation_out] = RetinalImageAlign_core(ImA,ImB,BlurSigma,CropB,rot_angles)
% RetinalImageAlign is the core function of the Align_Retinal_Images script
% (my testbed for this stuff)
%
% written by Niklas Domdei, April 2020

AllOpts = combvec(BlurSigma,CropB);
AllSteps = size(combvec(BlurSigma,CropB,rot_angles),2);
CStep = 1;

[BestSX,BestSY,BestError,BestRotation] = deal(zeros(size(AllOpts,2),1));
wbh = waitbar(0,'Testing alignment...');
for ax = 1:size(AllOpts,2)
    
    [ShiftX,ShiftY,Error] = deal(zeros(length(rot_angles),1));
    
    for bx = 1:length(rot_angles)

        waitbar(CStep/AllSteps,wbh)

        % crop to center
        CropVec = 1+(AllOpts(2,ax)):size(ImA,1)-(AllOpts(2,ax));
        ImA_cropped = ImA(CropVec,CropVec);
        ImA_cropped(ImA_cropped<0.05) = round(mean(ImA_cropped(:)));    % remove pure black (zeros) from image

        % rotate centered, then crop
        ImB_rotated = imrotate(ImB,rot_angles(bx),'crop');
        ImB_cropped_rotated = ImB_rotated(CropVec,CropVec);
        ImB_cropped_rotated(ImB_cropped_rotated<0.05) = round(mean(ImB_cropped_rotated(:)));

        %     figure,
        %     subplot(1,2,1)
        %     imshow(ImA_cropped)
        %     subplot(1,2,2)
        %     imshow(ImB_cropped_rotated)

        FFT2ImA = fft2(ImA_cropped);
        FFT2ImB = fft2(ImB_cropped_rotated);

        if AllOpts(1,ax)>0      % BlurSigma

            % blur of inner sphere (=low pass filtering)
            H = fspecial('gaussian',size(FFT2ImA,1),AllOpts(1,ax));
            H = H./max(H(:));
            H_inv = (H.*-1)+1;  % remove only the high frequencies (center of matrix)

            FFT2ImA_blurred = ifftshift(fftshift(FFT2ImA).*H_inv);
            FFT2ImA_blurred(1,1) = FFT2ImA(1,1);    % keep this value, it is extremely important for fft to work
            FFT2ImB_blurred = ifftshift(fftshift(FFT2ImB).*H_inv);
            FFT2ImB_blurred(1,1) = FFT2ImB(1,1);    % keep this value, it is extremely important for fft to work

            RefImFFT = FFT2ImA_blurred;
            CurImFFT = FFT2ImB_blurred;

            %     ImA_blurred = abs(ifft2(FFT2ImA_blurred));
            %     ImA_blurred = ImA_blurred./max(ImA_blurred(:));
            %     figure,
            %     subplot(2,2,1),imshow(ImA_cropped), axis square, caxis([0 1])
            %     subplot(2,2,2),imagesc(abs(fftshift(FFT2ImA))), axis square, caxis([0 300])
            %     subplot(2,2,3),imagesc(abs(fftshift(FFT2ImA_blurred))), axis square, caxis([0 300])
            %     subplot(2,2,4),imshow(ImA_blurred), axis square, caxis([0 1])

        else
            RefImFFT = FFT2ImA;
            CurImFFT = FFT2ImB;
        end

        % simple solution
        [output,~] = dftregistration(RefImFFT,CurImFFT,10);

        ShiftX(bx) = output(4);
        ShiftY(bx) = output(3);
        Error(bx) = output(1);
        
        CStep = CStep+1;
    end
    
    % clear nonsense data from vector
    NonSense = abs(ShiftX)>80 | abs(ShiftY)>80;
    ShiftX(NonSense) = [];
    ShiftY(NonSense) = [];
    Error(NonSense) = [];   
    CleanRA = rot_angles;
    CleanRA(NonSense) = [];   
    
    % store data with minimal error
    if ~isempty(Error)
        [~,BestOffsetIDX] = min(Error);
        BestSX(ax) = ShiftX(BestOffsetIDX);
        BestSY(ax) = ShiftY(BestOffsetIDX);
        BestError(ax) = Error(BestOffsetIDX);
        BestRotation(ax) = CleanRA(BestOffsetIDX);
    end

end
close(wbh)



%% analyse and sort the collected data using something similar to MAD (median absolute deviance) (like in find_PRL.m)
x = BestSX;
y = BestSY;
r = BestRotation;
outlierfactor = .66;  % pixel
outlierx = logical((x>median(x)+outlierfactor) + (x<median(x)-outlierfactor));
outliery = logical((y>median(y)+outlierfactor) + (y<median(y)-outlierfactor));
outliers = logical(outlierx+outliery);
xo = x; yo = y; ro = r;
xo(outliers) = [];
yo(outliers) = [];
ro(outliers) = [];

GoodMatch = 0;
if length(xo)>9
    GoodMatch = 1;
    disp('Alignment succeeded :)')
    
else
    disp('Alignment failed :(')
end



%% report median of the clean data back to calling function or script
ShiftX_out = median(xo);
ShiftY_out = median(yo);
Rotation_out = median(ro);

end