function[Second_Struct] = Struct_AddSegment(MRS_Struct,Second_Struct)
% CWJ: Stopgap function for copying segmentation results to several voxels

Second_Struct.mask = MRS_Struct.mask;
Second_Struct.out = MRS_Struct.out;

Fn = fieldnames(MRS_Struct.Fit.Metabolite);
tissuefra=MRS_Struct.Fit.Metabolite.(Fn{1}).SignalAmplitude/MRS_Struct.Fit.Metabolite.(Fn{1}).SignalAmplitude_CSF;
for J=1:length(Fn)
    Second_Struct.Fit.Metabolite.(Fn{J}).SignalAmplitude_CSF = Second_Struct.Fit.Metabolite.(Fn{J}).SignalAmplitude / tissuefra;
end

end