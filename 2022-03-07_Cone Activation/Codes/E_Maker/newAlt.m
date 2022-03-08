function [offset alt]=newAlt(intensity, correct, alt, trial, E_orientations) %#ok<*INUSD>

% CFG = getappdata(getappdata(0,'hAomControl'),'CFG');
load('BAK1041R_2020_04_20_09_53_26_AOSLO.mat')

if strcmp(CFG.optotype,'E')
    offset=intensity;
    alt=E_orientations(trial);  %random rotations for all trials are previously defined in E_orientations [0,90,180,270] degrees
    
    
    %alt=90*floor(rand*4);  %random rotation about [0,90,180,270] degrees
    
else
    offset=intensity;
    sign=rand-0.5;
    if sign<0
        offset=-offset;
        alt=1;
    else
        alt=0;
    end
end