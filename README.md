# Cone activation

MATLAB pipeline that estimates, for every cone, how strongly it is activated by a stimulus in AOSLO stabilised-video experiments (University Clinic Bonn).

The code is in `Codes/` (a copy under `2022-03-07_Cone Activation/`). The documentation with figures is `ReadMe.docx`.

## How to run

1. Register each individual sumNorm (from one video) to the master sumNorm (from the PRL videos). The master sumNorm provides the cone locations and CDCs. Use the registration code in `Codes/Niklas_Registration/`: it handles both shift and rotation and was the most accurate of the options tested.
2. If you do not have the stimuli saved, run `MyStim` from `Codes/E_Maker/`. It currently works for the tumbling E only. Its output `Es.mat` is specific to one experimental session; copy it into the folder with that session's videos.
3. Run `Codes/Last Version/mainAuto.m`. You are asked for the videos and data matrices at the start; then it runs by itself. Expect 20 to 30 s of computing per 1-second video.

## What mainAuto calls

- `QualityCheck_July14`: checks that the registration to the master sumNorm worked and that the cross is found in all frames that contain it. The cross is found from the shifts of offline stabilisation (`FindPRL_July14`, which only works for gain = 0 videos).
- `FineStimMaker_July16`: pads the stimulus image to a fixed size (the maximum in the session) and builds the array the pipeline needs.
- `SVA_oneStrip_traceAD`: the core of the pipeline. From the stimulus, cone locations and cross, it computes per-cone activation and stores it in `ConeLoc_actLogs.mat`. This holds the summed Gaussian cone probability distribution per cone (Voronoi cell). Uncomment line 250 to get the probability map itself.

## Time course

`timeCourse.m` computes the time course of the activations. It is not finalised: there were not enough examples to test it. The related papers are in `References_TimeCourse/`.

Older versions of the scripts are in `Codes/Previous Versions/`.
