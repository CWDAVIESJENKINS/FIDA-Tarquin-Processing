function[RawData]=DWMRS_ReadDat_Base(IN)

file_Met = IN;
RawData=io_loadspec_twix(file_Met);
Ndir = RawData.p.alFree{4};
Ng = RawData.p.alFree{5};
RepSize = Ndir*Ng+1;

RawData.p.RawDataLocation = file_Met; % Add data locations to structure

G = cell2mat(RawData.p.alFree(22:21+Ng));
RawData.Diff_Dir = [0 ,repmat(G,[1,Ndir])];
RawData.Diff_G = [0,reshape(repmat(cell2mat(RawData.p.alFree(15:15+Ndir-1)), [Ng,1]),[Ng*Ndir,1])'];

end