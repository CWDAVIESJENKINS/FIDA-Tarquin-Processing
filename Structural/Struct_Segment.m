function MRS_struct = Struct_Segment(MRS_struct,OutLoc)
% Runs segmentation script if segmented images not present according to
% file convention of c1, c2 and c3 as prefixes on the anatomical image name
% for the GM, WM and CSF segmentations. If these files are present, they
% are loaded and used for the voxel segmentation

% First check if SPM12 is installed and on the search path
spmversion = fileparts(which('spm'));
if isempty(spmversion)
    addpath('/cubric/software/spm')
    warning('Adding spm to Matlab search path')
end


%numscans = MRS_struct.p.numscans;

% Set up SPM for batch processing
spm('defaults','fmri');
spm_jobman('initcfg');


% ii nscans
% kk nvox

% 1 - Take nifti from GannetCoRegister and segment it in SPM
[T1dir, T1name, T1ext] = fileparts(MRS_struct.T1Reg.mask.T1image);
anatimage = MRS_struct.T1Reg.mask.T1image;

% Check to see if segmentation already done - if not, do it
tmp = [T1dir '/c1' T1name T1ext];
if ~exist(tmp,'file')
    disp('Starting Segmentation - SPM12')
    CallSPM12seg_Tarq(anatimage);
    disp('Segmentation Complete!')
end
% 2 - Determine GM, WM and CSF fractions for each voxel
keyboard
if strcmp(T1dir,'')
    T1dir = '.';
end

GM  = [T1dir '/c1' T1name T1ext];
WM  = [T1dir '/c2' T1name T1ext];
CSF = [T1dir '/c3' T1name T1ext];

GMvol  = spm_vol(GM);
WMvol  = spm_vol(WM);
CSFvol = spm_vol(CSF);

voxmaskvol = spm_vol(MRS_struct.T1Reg.mask.outfile);
[a,b,c] = fileparts(voxmaskvol.fname);

% GM
O_GMvox.fname = fullfile(a, [b '_GM' c]);
O_GMvox.descrip = 'GMmasked_MRS_Voxel_Mask';
O_GMvox.dim = voxmaskvol.dim;
O_GMvox.dt = voxmaskvol.dt;
O_GMvox.mat = voxmaskvol.mat;
GM_voxmask_vol = GMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
O_GMvox = spm_write_vol(O_GMvox, GM_voxmask_vol);

% WM
O_WMvox.fname = fullfile(a, [b '_WM' c]);
O_WMvox.descrip = 'WMmasked_MRS_Voxel_Mask';
O_WMvox.dim = voxmaskvol.dim;
O_WMvox.dt = voxmaskvol.dt;
O_WMvox.mat = voxmaskvol.mat;
WM_voxmask_vol = WMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
O_WMvox = spm_write_vol(O_WMvox, WM_voxmask_vol);

% CSF
O_CSFvox.fname = fullfile(a, [b '_CSF' c]);
O_CSFvox.descrip = 'CSFmasked_MRS_Voxel_Mask';
O_CSFvox.dim = voxmaskvol.dim;
O_CSFvox.dt = voxmaskvol.dt;
O_CSFvox.mat = voxmaskvol.mat;
CSF_voxmask_vol = CSFvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
O_CSFvox = spm_write_vol(O_CSFvox, CSF_voxmask_vol);

%CWJ Whole brain
O_WB.fname = fullfile(a, [b '_WholeBrain' c]);
O_WB.descrip = 'WholeBrain_Mask';
O_WB.dim = voxmaskvol.dim;
O_WB.dt = voxmaskvol.dt;
O_WB.mat = voxmaskvol.mat;
WB_voxmask_vol = logical(GMvol.private.dat(:,:,:) + WMvol.private.dat(:,:,:) + CSFvol.private.dat(:,:,:));
O_WB = spm_write_vol(O_WB, WB_voxmask_vol);

% 3 - Calculate an adjusted iu and output it to the structure

gmsum = sum(sum(sum(O_GMvox.private.dat(:,:,:))));
wmsum = sum(sum(sum(O_WMvox.private.dat(:,:,:))));
csfsum = sum(sum(sum(O_CSFvox.private.dat(:,:,:))));

gmfra = gmsum/(gmsum+wmsum+csfsum);
wmfra = wmsum/(gmsum+wmsum+csfsum);
csffra = csfsum/(gmsum+wmsum+csfsum);

%tissuefra = gmfra+wmfra;

MRS_struct.T1Reg.tissue.GMfra = gmfra;
MRS_struct.T1Reg.tissue.WMfra = wmfra;
MRS_struct.T1Reg.tissue.CSFfra = csffra;

% 4 - Build output
if ishandle(104)
    clf(104); % MM (170831)
end
h = figure(104);
% MM (170831): Open figure in center of screen
scr_sz = get(0, 'ScreenSize');
fig_w = 1000;
fig_h = 707;
set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
set(h,'Color',[1 1 1]);
figTitle = 'Segment Output';
set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');

% Output results
subplot(2,3,4:6);
axis off;

text_pos = 0.87;

tmp1 = 'GM voxel fraction: ';
tmp2 = sprintf(' %.3g', MRS_struct.T1Reg.tissue.GMfra);
text(0.5, text_pos-0.12, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
text(0.5, text_pos-0.12, tmp2, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13);

tmp1 = 'WM voxel fraction: ';
tmp2 = sprintf(' %.3g', MRS_struct.T1Reg.tissue.WMfra);
text(0.5, text_pos-0.24, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
text(0.5, text_pos-0.24, tmp2, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13);

tmp1 = 'CSF voxel fraction: ';
tmp2 = sprintf(' %.3g', MRS_struct.T1Reg.tissue.CSFfra);
text(0.5, text_pos-0.36, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
text(0.5, text_pos-0.36, tmp2, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13);

tmp1 = 'Filename: ';
[~,tmp2,tmp3] = fileparts(MRS_struct.p.RawDataLocation);
text(0.5, text_pos-0.48, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
text(0.5, text_pos-0.48, [' ' tmp2 tmp3],  'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');

tmp1 = 'Anatomical image: ';
[~,tmp2,tmp3] = fileparts(MRS_struct.T1Reg.mask.T1image);
text(0.5, text_pos-0.6, tmp1, 'FontName', 'Helvetica', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
text(0.5, text_pos-0.6, [' ' tmp2 tmp3],  'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');


% Voxel segmentation (MM: 180807)
T1 = spm_read_vols(spm_vol(anatimage));
img_t     = flipud(voxel2world_space(spm_vol(anatimage), MRS_struct.p.voxoff));
vox_t     = flipud(voxel2world_space(voxmaskvol, MRS_struct.p.voxoff));
vox_t_GM  = flipud(voxel2world_space(O_GMvox, MRS_struct.p.voxoff));
vox_t_WM  = flipud(voxel2world_space(O_WMvox, MRS_struct.p.voxoff));
vox_t_CSF = flipud(voxel2world_space(O_CSFvox, MRS_struct.p.voxoff));
img_t = img_t/max(T1(:));
img_montage = [img_t+0.175*vox_t, img_t+0.21*vox_t_GM, img_t+0.25*vox_t_WM, img_t+0.4*vox_t_CSF];

h = subplot(2,3,1:3);
imagesc(img_montage);
colormap('gray');
img = MRS_struct.T1Reg.mask.img;
img = img(:);
caxis([0 mean(img(img>0.01)) + 3*std(img(img>0.01))]);
axis equal;
axis tight;
axis off;
text(floor(size(vox_t,2)/2), 20, 'Voxel', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center'); 
text(floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 20, 'GM', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center'); 
text(2*floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 20, 'WM', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center'); 
text(3*floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 20, 'CSF', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center'); 
get(h,'pos');
set(h,'pos',[0 0.15 1 1]);


[~,tmp,tmp2] = fileparts(MRS_struct.p.RawDataLocation);
[~,tmp3,tmp4] = fileparts(MRS_struct.T1Reg.mask.T1image);
t = ['Voxel from ' tmp tmp2 ' on ' tmp3 tmp4];
title(t, 'FontName', 'Helvetica', 'FontSize', 15, 'Interpreter', 'none');

set(gcf,'PaperUnits','inches');set(gcf,'PaperSize',[11 8.5]);set(gcf,'PaperPosition',[0 0 11 8.5]);

[path,name] = fileparts(MRS_struct.p.RawDataLocation);
if ~exist('OutLoc','var')
    CJ_SavFig(gcf,[name,'_Seg'],path);
else
    CJ_SavFig(gcf,[name,'_Seg'],OutLoc);
end

end

