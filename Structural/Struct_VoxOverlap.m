function MRS_struct = Struct_VoxOverlap(Struct1 , Struct2 , Nii1, Nii2)
% function MRS_struct = Struct_CoRegister(MRS_struct, nii_name , OutLoc)
%
% CWJ: Function for the coregistration of MRS voxels using SPM calls. Based
% on Gannet code for GM/WM segmentation. Function extracts voxel position 
% using parameters in twix data. A  mask nii is generated at the specified 
% location
% Input    MRS_struct = Structure currently used by my raw data readers
%          nii_name = Full path to nii files
% Optional OutLoc = Directory to save output nii and figure

% First check if SPM12 is installed and on the search path
spmversion = fileparts(which('spm'));
if isempty(spmversion)
    addpath('/cubric/software/spm')
    warning('Adding spm to Matlab search path')
end
warning('off','MATLAB:nearlySingularMatrix');warning('off','MATLAB:qhullmx:InternalWarning');

S1 = CoReg(Struct1, Nii1);
S2 = CoReg(Struct2, Nii2);

keyboard


end

%%
function[MRS_struct] = CoReg(MRS_struct , nii_name , OutLoc)
% Define location of ouput mask, if none specified
[path,name] = fileparts(MRS_struct.p.RawDataLocation);
if ~exist('OutLoc','var')
    fidoutmask = fullfile(path,[name '_mask.nii']);
else
    fidoutmask = fullfile(OutLoc,[name '_mask.nii']);
end

% Extract voxel position and rotation parameters from MRS_struct (from FID-A structure)
NormSag = MRS_struct.p.NormSag;NormCor = MRS_struct.p.NormCor;NormTra = MRS_struct.p.NormTra;
VoI_InPlaneRot = MRS_struct.p.VoI_InPlaneRot;

% Correct voxel offsets by table position (if field exists)
if isfield(MRS_struct.p,'TablePosition') && ~isempty(MRS_struct.p.TablePosition)
    VoxOffs = [MRS_struct.p.voxoff(1)+MRS_struct.p.TablePosition(1) MRS_struct.p.voxoff(2)+MRS_struct.p.TablePosition(2) MRS_struct.p.voxoff(3)+MRS_struct.p.TablePosition(3)];
else
    VoxOffs = [MRS_struct.p.voxoff(1) MRS_struct.p.voxoff(2) MRS_struct.p.voxoff(3)];
end


% Parse direction cosines of the MRS voxel's normal vector and the rotation angle
% around the normal vector
% The direction cosine is the cosine of the angle between the normal
% vector and the respective direction.
% Example: If the normal vector points exactly along the FH direction, then: 
% NormSag = cos(90) = 0, NormCor = cos(90) = 0, NormTra = cos(0) = 1.
Norm = [-NormSag -NormCor NormTra];
ROT = VoI_InPlaneRot;
% Find largest element of normal vector of the voxel to determine primary
% orientation. 
% Example: if NormTra has the smallest out of the three Norm
% values, the angle of the normal vector with the Tra direction (FH) is the
% smallest, and the primary orientation is transversal.
[~, maxdir] = max([abs(NormSag) abs(NormCor) abs(NormTra)]);
switch maxdir
    case 1
        vox_orient = 's'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
    case 2
        vox_orient = 'c'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
    case 3
        vox_orient = 't'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
end
    
% Phase reference vector
% Adapted from Rudolph Pienaar's "vox2ras_rsolveAA.m" and
% Andre van der Kouwe's "autoaligncorrect.cpp"
Phase	= zeros(3, 1);
switch vox_orient
    case 't'
        % For transversal voxel orientation, the phase reference vector lies in
        % the sagittal plane
        Phase(1)	= 0;
        Phase(2)	=  Norm(3)*sqrt(1/(Norm(2)*Norm(2)+Norm(3)*Norm(3)));
        Phase(3)	= -Norm(2)*sqrt(1/(Norm(2)*Norm(2)+Norm(3)*Norm(3)));
        VoxDims = [MRS_struct.p.voxdim(1) MRS_struct.p.voxdim(2) MRS_struct.p.voxdim(3)];
    case 'c'
        % For coronal voxel orientation, the phase reference vector lies in
        % the transversal plane
        Phase(1)	=  Norm(2)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(2)	= -Norm(1)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(3)	= 0;
        VoxDims = [MRS_struct.p.voxdim(1) MRS_struct.p.voxdim(2) MRS_struct.p.voxdim(3)];
    case 's'
        % For sagittal voxel orientation, the phase reference vector lies in
        % the transversal plane
        Phase(1)	= -Norm(2)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(2)	=  Norm(1)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(3)	= 0;
        VoxDims = [MRS_struct.p.voxdim(1) MRS_struct.p.voxdim(2) MRS_struct.p.voxdim(3)];
end

% The readout reference vector is the cross product of Norm and Phase
Readout = cross(Norm, Phase);
M_R = zeros(4, 4);
M_R(1:3, 1)	= Phase;
M_R(1:3, 2)	= Readout;
M_R(1:3, 3) = Norm;

% Define matrix for rotation around in-plane rotation angle
M3_Mu	= [	 cos(ROT)	sin(ROT)	0
            -sin(ROT)	cos(ROT)	0
            0           0           1];
        
M3_R	= M_R(1:3,1:3)	* M3_Mu;
M_R(1:3,1:3)	= M3_R;

% The MGH vox2ras matrix inverts the Readout column
M_R		= M_R *   [ 1  0  0  0
                    0 -1  0  0
                    0  0  1  0
                    0  0  0  1];

% Final rotation matrix
rotmat = M_R(1:3,1:3);

V = spm_vol(nii_name);
[T1,XYZ] = spm_read_vols(V);

%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim] = spm_get_bbox(V,'fv'); % MM (180220)
voxdim = abs(voxdim)';
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

% We need to flip ap and lr axes to match NIFTI convention
VoxOffs(1) = -VoxOffs(1);
VoxOffs(2) = -VoxOffs(2);

% Define voxel coordinates before rotation and transition
vox_ctr = ...
    [VoxDims(1)/2 -VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2 -VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2  VoxDims(2)/2  VoxDims(3)/2;
     VoxDims(1)/2  VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2  VoxDims(2)/2 -VoxDims(3)/2;
     VoxDims(1)/2  VoxDims(2)/2 -VoxDims(3)/2;
     VoxDims(1)/2 -VoxDims(2)/2 -VoxDims(3)/2;
    -VoxDims(1)/2 -VoxDims(2)/2 -VoxDims(3)/2];

% Apply rotation as prescribed
vox_rot = rotmat*vox_ctr.';

% Shift rotated voxel by the center offset to its final position
vox_ctr_coor = [VoxOffs(1) VoxOffs(2) VoxOffs(3)];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot + vox_ctr_coor;

% Create a mask with all voxels that are inside the VOI
mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((VoxDims(1)/2)^2+(VoxDims(2)/2)^2+(VoxDims(3)/2)^2);
distance2voxctr = sqrt(sum((XYZ-repmat([VoxOffs(1) VoxOffs(2) VoxOffs(3)].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr <= sphere_radius) = 1;
mask(sphere_mask == 1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);
tri = delaunayn([vox_corner.'; [VoxOffs(1) VoxOffs(2) VoxOffs(3)]]);
tn = tsearchn([vox_corner.'; [VoxOffs(1) VoxOffs(2) VoxOffs(3)]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

% Take over the voxel dimensions from the structural
mask = reshape(mask, V.dim);

V_mask.fname   = fidoutmask ;
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim     = V.dim;
V_mask.dt      = V.dt;
V_mask.mat     = V.mat;

V_mask = spm_write_vol(V_mask,mask);

% Build output page
%fidoutmask = cellstr(fidoutmask);
MRS_struct.T1Reg.mask.outfile(1,:) = fidoutmask;
% Not clear how to formulate the rotations for triple rotations (revisit)
%MRS_struct.p.voxang(1,:) = [NaN NaN NaN];  

% Transform structural image and co-registered voxel mask from voxel to
% world space for output (MM: 180221)
[img_t,img_c,img_s] = voxel2world_space(V,VoxOffs);
[mask_t,mask_c,mask_s] = voxel2world_space(V_mask,VoxOffs);

img_t = flipud(img_t/max(T1(:)));
img_c = flipud(img_c/max(T1(:)));
img_s = flipud(img_s/max(T1(:)));

img_t = img_t + 0.175*flipud(mask_t);
img_c = img_c + 0.175*flipud(mask_c);
img_s = img_s + 0.175*flipud(mask_s);

size_max = max([max(size(img_t)) max(size(img_c)) max(size(img_s))]);
three_plane_img = zeros([size_max 3*size_max]);
three_plane_img(:,1:size_max)              = image_center(img_t, size_max);
three_plane_img(:,size_max+(1:size_max))   = image_center(img_s, size_max);
three_plane_img(:,size_max*2+(1:size_max)) = image_center(img_c, size_max);

MRS_struct.T1Reg.mask.img = three_plane_img;
MRS_struct.T1Reg.mask.T1image(1,:) = nii_name;

warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');
keyboard
end
