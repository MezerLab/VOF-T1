function vofGetCandidates(wholebrainfgPath,fsROIdir,outdir,thresh,v_crit)
% Segment the VOF from a wholebrain connectome
%
% [L_VOF, R_VOF, L_pArc, R_pArc, L_pArc_vot, R_pArc_vot] = AFQ_FindVOF(wholebrainfgPath,L_arcuate,R_arcuate,fsROIdir,outdir,thresh,v_crit, dt)
%
% This function is a modified version of AFQ's function AFQ_FindVOF.m
% (Copyright Jason D. Yeatman, September 2014. Code released with: Yeatman
% J.D., Weiner K.S., Pestilli F., Rokem A., Mezer A., Wandell B.A. (2014).
% The vertical occipital fasciculus: A forgotten highway. PNAS)
%
% It is used to extract and save the vertical candidate streamlines
% (Left_VerticalFG.mat and Right_VerticalFG.mat).
%
% Inputs:
%
% wholebrainfgPath - A path (or fg structure) for a wholebrain fiber group.
% L_arcuate        - Segmented arcuate fasciculus (left hemisphere). See 
%                    AFQ_SegmentFiberGroups
% R_arcuate        - Segmented arcuate fasciculus (left hemisphere).
% fsROIdir         - Path to a directory containing .mat ROIs of each
%                    cortical region that is segmnted by freesurfer. This
%                    means that you must first run freesurfers recon-all on
%                    a t1 weighted image to get a cortical segmentation.
%                    Next use the function: 
%                    fs_roisFromAllLabels(fsIn,outDir,type,refT1)
%                    to convert the freesurfer segmentation into ,mat ROIs
% outdir           - This is where all the outputs will be saved
% thresh           - A fiber must travel vertical for a large proportion of
%                    its length. The default values are likely fine
% vcrit            - To find fibers that we can considder vertical, we must
%                    define a threshold of how far a fiber must travel in
%                    the vertical direction. v_crit defines how much
%                    farther a fiber must travel in the vertical direction
%                    compare to other directions (e.g., longitudinal) to be
%                    a candidate for the VOF. The default is 1.3
%
% Outputs
% Left_VerticalFG.mat  - Left  hemisphere candidate streamlines fiber groups
% Right_VerticalFG.mat - Right hemisphere candidate streamlines fiber groups

%% Argument checking and parameter setting
% Path to ROIs
if notDefined('fsROIdir')
    fsROIdir = uigetdir([],'Select an ROI directory');
end
% output directory
if notDefined('outdir')
    outdir = uigetdir([],'Select an output directory');
end
% Remove any fiber that doesn't go vertical (positive z) for thresh% of its
% coordinates
if notDefined('thresh')
    thresh = [.95 .6];
end
% Fibers must travel this much farther vertically than other directions
if notDefined('v_crit')
    v_crit = 1.3;
end

%% Find vertical fibers
% From the wholebrain fiber group find all the vertical fibers that
% terminate in ventral occipitotemporal cortex (VOT).
[L_fg_vert, R_fg_vert] = AFQ_FindVerticalFibers(wholebrainfgPath,fsROIdir,outdir,thresh,v_crit);

%% Remove streamlines shorter than 20 mm
% Note: this usually happens automatically in AFQ_FindVerticalFibers,
% however if the resolution is smaller than 2 mm, and specifically if the
% step size is smaller than 1 mm.
stepSize = median(vecnorm(L_fg_vert.fibers{2}(:,1:end-1)-L_fg_vert.fibers{2}(:,2:end)));
% Left
fibLen = cellfun(@length,L_fg_vert.fibers);
minLen = ceil(20/stepSize);
indicesKeep = find(fibLen>=minLen);
fg = fgRetainIndices(L_fg_vert,indicesKeep);
% Right
fibLen = cellfun(@length,R_fg_vert.fibers);
minLen = ceil(20/stepSize);
indicesKeep = find(fibLen>=minLen);
fg = fgRetainIndices(R_fg_vert,indicesKeep);


%% Save to files
dtiWriteFiberGroup(L_fg_vert,fullfile(outdir,'Left_VerticalFG'));
dtiWriteFiberGroup(R_fg_vert,fullfile(outdir,'Right_VerticalFG'));
