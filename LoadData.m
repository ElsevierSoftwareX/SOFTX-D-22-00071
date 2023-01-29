%#################################################################################################################
%PURPOSE :Extracts data and times from data files.
%VERSION: 1.0 [10/04/2017]
%         1.1 [20/07/2020] (vb) Application of a filter to erase initial distorted data points.
%         1.2 [31/05/2021] (vb) Added the explicit filepath.
%DATE    :10/04/2017
%CHANGES :1.0 [] 
%AUTHOR  :VB.
%##################################################################################################################
function [CommentTS, N_T1, N_T2, t_T1, t_T2, S] = LoadData(Data_Folder, DataFileName, TimeRowFileName, TimeColumnFileName, EraseCol, EraseRow)
 %
 TimeRowFileName=[Data_Folder TimeRowFileName];
 fid = fopen(TimeRowFileName);  %
  t_T1 = fscanf(fid,'%f');
 fclose(fid);
 %N_T1=size(t_T1);
 %
 TimeColumnFileName=[Data_Folder TimeColumnFileName];
 fid = fopen(TimeColumnFileName);  %
  t_T2 = fscanf(fid,'%f');
 fclose(fid);
 %N_T2=size(t_T2);
 %
 DataFileName=[Data_Folder DataFileName];
 S = dlmread(DataFileName);
 %
 %Erasing Filter
 if (EraseCol||EraseRow)~=0
   S(:, 1:EraseCol)=[];
   S(1:EraseRow, :)=[];
   t_T1(1:EraseRow)=[];
   t_T2(1:EraseCol)=[];
 end
 %
 CommentTS='';
 N_T1=size(t_T1);
 N_T2=size(t_T2);
 return;
%
end

