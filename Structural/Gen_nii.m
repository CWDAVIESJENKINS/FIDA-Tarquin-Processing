function [ Loc , CMDout ] = Gen_nii( Dcm , Nii  )
% Function for generating nifti files using dicom data. Uses CUBRIC dcm2nii
% script, and saves .nii in dicom directory subfolder "*/nii"
% Input:    Target = Target directory or dicom
% Optional: Nii = Location to save Nii
% Output:   Loc = location of first nifti file in target directory
%           CMDout = Command line output from dcm2nii (debugging)

[Path,~,~] = fileparts(Dcm); % If a directory is not input, assume parent directory is target

if ~exist('Nii','var')
    Nii = [Path , '/nii']; % Delineate target directory for .nii files if none specified
    if ~exist(Nii,'dir')
        mkdir(Nii)
    end
end

RUN = ['dcm2niix ', Nii, ' ', Path];
[ ~ , CMDout] = system(RUN); % Run cubric script to extract .nii

List = dir([Nii,'/*.gz']);

for J=1:length(List)
    gunzip([Nii,'/',List(J).name], Nii); delete([Nii,'/',List(J).name]);
end

List = dir([Nii,'/*.nii']); Loc = [List(1).folder , '/' , List(1).name]; % State locaton of first .nii for use in CoRegister


end