function[] = DWMRS_NonLinearity(Specs,Mask_Nii,Dir)


%% gradnonlin correction
scanner = 'connectom';

obj.image_file = '/cubric/data/sapcj3/dMRS/20200303/meas_MID394_dMRS_gmax300_5grad_3dir_met_FID10665_mask.nii';
obj.mask_file = obj.image_file;

[ ~, obj.dev_file ] = func_create_disco_file( obj.image_file, scanner );
[ ppd, nnd, extd ]  = Supp_gz_fileparts( obj.image_file );
def_file = [ ppd filesep nnd '__deform_grad_rel' extd ];

%% compute the voxel-specific b-matricies

obj.eddypars = '/cubric/data/sapcj3/dMRS/20200303/220211-1.eddy_parameters'; % this file should have 16 columns and number of bvalues * number of directions columns
obj.bvec_file_init = '/cubric/data/sapcj3/dMRS/20200303/20200303.bvec'; % the directions, now set to [0 0 1], [-1 0 0] and [0 1 0] is that still correct?
obj.bval_file = '/cubric/data/sapcj3/dMRS/20200303/20200303.bval'; % the new b-values that Elena provided
par = []; % can remain empty
[ obj.mod_x, obj.mod_y, obj.mod_z ] = DWI_preproc_disco_grad_dev( obj, par );
[modx,~] = mdm_nii_read(obj.mod_x); [mody,~] = mdm_nii_read(obj.mod_y); [modz,hz] = mdm_nii_read(obj.mod_z);
bc = sqrt(modx.^2+mody.^2+modz.^2);
bval_c_file = '/cubric/data/sapcj3/dMRS/20200303/corrected_bvalues.nii.gz'; %adapt this

mdm_nii_write(bc,bval_c_file,hz);

%% plot imposed vs corrected b-values
% mask_DISCO_file = '/cubric/data/sapct5/DW_MRS/mask_s_DISCO.nii.gz';
mask_DISCO_file = obj.mask_file;
[msk,h] = mdm_nii_read(mask_DISCO_file);
b = load(obj.bval_file);
msk(msk<0) = 0; msk = single(msk);

% each voxel in mask contributes equally
bc1_m = mean(vec(bc,logical(msk)),2);
figure; scatter(b,bc1_m,'b'); xlabel('imposed b-value'); ylabel('corrected b-value'); hold all
plot([min(b);max(b)],[min(b);max(b)],'--k'); axis square

% each voxel in mask contributes according to fraction
%bc2_m = repmat(msk,[1 1 1 size(bc,4)]).*bc;
%bc2_m = mean(vec(bc2_m,logical(msk)),2);
%scatter(b,bc2_m,'r'); 

%dlmwrite('/cubric/data/sapct5/DW_MRS/1408/corrected_bvalues_mask.bval',bc1_m')