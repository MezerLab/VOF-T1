function fg = fgFindTerminatingAtCoretx(fg,aparcasegFile)
% fgFindTerminatingAtCoretx takes as input a fiber group
% fg and a aparc+aseg.nii.gz FreeSurfer file, and finds streamlines that 
% terminate within 2 mm from the cortex (it could be modified to
% include subcortical gray matter as well).
%
% This is useful for cases in which the initial tractography allows
% streamlines to terminate far from the cortex. For example, using MRtrix's
% ACT tractography sometimes includes streamlines terminating near the
% ventricles, as these regions are considered as GM/WM interface by the 5tt
% algorithm.
%
% The indices of streamlines that terminate near the cortex are saved in a
% param field in the fg structure, under the name "terminating_prematurely".

tmp = readFileNifti(aparcasegFile);
tmp.data(tmp.data<1000) = 0;
cortexFile = fullfile(fileparts(aparcasegFile),'aparc+aseg_gm.nii.gz');
dtiWriteNiftiWrapper(double(logical(tmp.data)),tmp.qto_xyz,cortexFile);


fgTmp = fg;
fgTmp.fibers = cellfun(@(x) x(:,end), fgTmp.fibers,'UniformOutput',false); % only save first and last nodes
targetROIfile = cortexFile;
targetROI2file = targetROIfile;
thresholdmm = 2;
[~, keep] = dtiSegmentFiberWithNiftiRoi(fgTmp, targetROIfile, targetROI2file, thresholdmm); % This could take a few minutes

% Add a params field
fg = fgRemoveParams(fg,'not_terminating_prematurely');
perPointFlag = 0;
fg.params{end+1}.name = 'not_terminating_prematurely';
fg.params{end}.stat = keep;