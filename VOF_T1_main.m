% This is the example code and data associated with the following paper:
% Schurr, Roey, Shir Filo, and Aviv A. Mezer. "Tractography delineation of the vertical occipital fasciculus using quantitative T1 mapping." NeuroImage (2019): 116121.‏

%% Set data dir
exampleDataDir = '/ems/elsc-labs/mezer-a/code/roey/vof_github/ExampleData/WH008';

%% Preprocessing
%% (1.1) Run FreeSurfer
% Use FreeSurfer's function recon-all to an anatomical T1-weighted image.
% Preferably, it should be aligned with the diffusion data. Otherwise, the
% out put will have to be registered to the diffusion data later.

%% (1.2) Convert the FreeSurfer ROIs to VistaSoft's .mat format
% Use fs_roisFromAllLabels(fsIn,outDir,type,refT1)

%% (1.3) (optional) Extract the MFS label
% If using, please cite Weiner et al. (2018): "On object selectivity and
% the anatomy of the human fusiform gyrus"

%% (2.1) Run whole-brain tractography
% Consider using MRtrix3's probabilistic tractogrpahy with ACT.


%% Post-processing
% Choose hemisphere
hemi = 'right'; % 'right' or 'left'
if strcmp(hemi,'left')
    azEl = [-60 30]; %azimuth and elavation for plot
else
    azEl = [60 30];
end

%% (1) Get vertical candidate streamlines
% The vertical candidates will be saved under outdir (Left_VerticalFG.mat and Right_VerticalFG.mat).
wholebrainfgFile = fullfile(exampleDataDir,'Diffusion','tck500000_act_backtrack_crop_at_gmwmi.mat');

fsROIdir = fullfile(exampleDataDir,'FreeSurfer/ROIs');
outdir = fullfile(exampleDataDir,'VOF');
if ~exist(outdir,'dir')
    mkdir(outdir)
end
thresh = [.95 .6];
v_crit = 1.3;
vofGetCandidates(wholebrainfgFile,fsROIdir,outdir,thresh,v_crit)

if strcmp(hemi,'right')
    fg = fgRead(fullfile(outdir,'Right_VerticalFG.mat'));
else
    fg = fgRead(fullfile(outdir,'Left_VerticalFG.mat'));
end

%% (2) (optional) Eliminate streamlines terminating prematurely
% The indices of streamlines that terminate prematurely are saved in a param
% field in the fg structure, under the name "terminating_prematurely".
% This is important if the seeding mask (in this case a WM/GM-interface
% mask from MRtrix3's 5ttgen) includes the walls of the ventricles.
aparcasegFile = fullfile(exampleDataDir,'FreeSurfer','aparc+aseg.nii.gz');
fg = fgFindTerminatingAtCoretx(fg,aparcasegFile);
keep = find(fgGetParams(fg,'not_terminating_prematurely'));
fg = fgRetainIndices(fg,keep);

%% (3) (optional) Ignore streamlines with too much extrapolated B1+
ignoreB1extrapolated = 1;
if ignoreB1extrapolated
    B1File = fullfile(exampleDataDir,'T1','B1_noExtrapolaiton_2DTI_mask.nii.gz'); % A mask of voxels in which the B1+ field was estimated directly, and not extrapolated
    wmFile = fullfile(exampleDataDir,'T1','wm_prob_2DTI.nii.gz'); % A white-matter probability map created using MRTrix's function 5ttgen
    
    wm = readFileNifti(wmFile);
    B1_img = readFileNifti(B1File);
    B1_img.data(wm.data<1) = nan; % Mask B1 to WM
    perPointFlag = 1;
    fg = dtiCreateQuenchStats(fg, 'B1mask', 'B1mask', perPointFlag, B1_img,'nanmedian'); % this is my version of dtiCreateQuenchStats, where in line 214 a "std" options exists.
    
    B1profiles = {fg.pathwayInfo(:).point_stat_array};
    B1profiles = cellfun(@(x) x(end,:), B1profiles, 'UniformOutput',false);
    noB1prcnt = cellfun(@(x) length(find(isnan(x)))/length(x), B1profiles); % The percent of nodes in each streamline, in which B1 was extrapolated (and not estimated or interpolated)
    
    fg.params{end+1}.name = 'noB1prcnt';
    fg.params{end}.stat = noB1prcnt;
end
%% (4) Extract VOF based on its T1-Mdn profile
if strcmp(hemi,'right')
    mfsFile = fullfile(exampleDataDir,'FreeSurfer','rh_MFS_2DTI.nii.gz');
else
    mfsFile = fullfile(exampleDataDir,'FreeSurfer','lh_MFS_2DTI.nii.gz');
end
paramImgFile = fullfile(exampleDataDir,'T1','T1_map_Wlin_2DTI.nii.gz');
paramImg = readFileNifti(paramImgFile);
% Ignore voxels outside the white matter mask
wmFile = fullfile(exampleDataDir,'T1','wm_prob_2DTI.nii.gz'); % A white-matter probability map created using MRTrix's function 5ttgen
wm = readFileNifti(wmFile);
paramImg.data(wm.data<1) = nan;
% Perform the separation
vofIndices = vofSeparateFromCandidates(fg,paramImg,mfsFile);
fgVof = fgRetainIndices(fg,vofIndices);

%% (5) Plot candidate streamlines colored by their T1-Mdn
r1 = readFileNifti(paramImgFile);
r1.data = 1./r1.data;
r1.data(isinf(r1.data)) = 0;
r1.data(r1.data<0.25) = 0;
r1.data(r1.data>1.25) = 0;

matplotlib_colormaps; % Load the color maps
colormap = plasmadata; % magmadata plasmadata
perPointFlag = 0;

fg = dtiCreateQuenchStats(fg, 'param_median', 'param', perPointFlag, paramImg, 'nanmedian', 1); % This is a version of Vistasoft's dtiCreateQuenchStats, where in line 218 a "nanstd" options exists
vals = fgGetParams(fg, 'param_median');
vals(vals>prctile(vals,90)) = prctile(vals,90); % This is what I usually use for T1-Mdn
vals(vals<prctile(vals,5)) = prctile(vals,5); % This is what I usually use for T1-Mdn

rgb = vals2colormap(vals,colormap);

AFQ_RenderFibers(fg, 'color', rgb,'tubes',[1],'radius',[0.3,1],'jittershading',0,'camera',azEl)
AFQ_AddImageTo3dPlot(r1, [0,0,-10]);
AFQ_AddImageTo3dPlot(r1, [1,0,0]);

%% (6) Plot the final VOF
AFQ_RenderFibers(fgVof, 'color', lines(1),'tubes',[1],'radius',[0.3,1],'jittershading',0.3,'camera',azEl)
AFQ_AddImageTo3dPlot(r1, [0,0,-10]);
AFQ_AddImageTo3dPlot(r1, [1,0,0]);