%########################################################################################################################################
%NAME    : main.m
%PURPOSE : Main program interface. 
%
%NOTES   : for T1-T2 data it is used a IR-CPMG kernel or SR-CPMG, for T2-T2 data it is used a 
%          CPMG-CPMG  kernel.
%#########################################################################################################################################
clear all;
close all; 
clc;
parFit=struct;
parFile=struct;
parFL=struct;
fprintf('MUpen2D Starts \r\n');
%
selpath = uigetdir('./DATA/','Open Data Directory');
%
endout=regexp(selpath,filesep,'split');
FolderName=endout{end};
Test_folder=[FolderName '/'];
Data_folder=['./DATA/' Test_folder];%addpath(Data_Folder);
Out_folder='./output_files/';%addpath(Out_folder);
 %-------------------------------------------------------------------------
NameFileFlags=[Data_folder 'FileFlag.par'];
NameFileSetInput=[Data_folder 'FileSetInput.par'];
NameFilePar=[Data_folder 'FilePar.par'];
%
[parFile]=SetInputFile(NameFileSetInput, parFile, 0);
%############################### FLag parameters #######################################
% FL_typeKernel=1;          %1 IR-CPMG; 4 T2-T2
% FL_InversionTimeLimits=0; %1 automatic, 0 manually selection inversion time ranges 
% FL_OutputData=0;          %1 create output data file for ILT2D
%
% load flag fom file
[parFL]= LoadFlags(NameFileFlags,0);  
%
%#################################### Load Data Set ####################################
%(vb)[31/05/2021] Added the explict filepath.
[CommentTS, N_T1, N_T2, Tau1, Tau2, S] = LoadData(Data_folder, parFile.filenamedata, parFile.filenameTimeY,parFile.filenameTimeX, parFL.FL_EraseCol, parFL.FL_EraseRow);
% 
B=eye(size(S,2)); %now not used, for future implementation
%############################# Set the Parameter structure ##########################
[parFit]=SetPar(NameFilePar,S, B, parFit,0);
parFit.Verbose=parFL.FL_Verbose;
parFit.Debug=parFL.FL_Debug;
parFit.Kernel=parFL.FL_typeKernel;
parFit.Amp_scale=parFL.FL_Amp_scale;
parFit.scale_fact= parFL.FL_Scale_fact;
%
h=msgbox('Please wait, computation can take a while ...','MUpen2DTool is running','warn');
try 
 % computation here % 
  [x,T1,T2,hist]=MUpen2D(Test_folder, Data_folder, Out_folder, parFit, parFile, parFL, CommentTS, N_T1, N_T2, Tau1, Tau2, S);
  if(parFL.FL_Offset==1)
   x(:, end)=[];
   x(end,:)=[];
  end
  grafico_1D(x,T1,T2, '1D Distribution', parFL.FL_typeKernel, Out_folder);
  grafico_2D(x,T1,T2,0 ,parFL.FL_typeKernel, '2D Map');
  grafico_2D(x,T1,T2,parFL.FL_NoContour ,parFL.FL_typeKernel, '2D Map Cotour');
  grafico_3D(x,T1,T2,parFL.FL_typeKernel, '3D Distribution');
  grafico_3D_data(hist.X,Tau1,Tau2,parFL.FL_typeKernel, 'Computed Data');
  grafico_3D_data(S,Tau1,Tau2,parFL.FL_typeKernel, 'Input Data');
  %------------------------------------------------------------------------
  % Data Analysis
  %
  if(parFL.FL_Stat)
    out_data=Residual_Analysis(hist.X, S, parFL);
      fid=fopen('./output_files/Parameters.txt','a');
      fprintf(fid,'\n===================================================================\n');
      fprintf(fid,'    Statistical Analysis \n');
      fprintf(fid,' Normal  Distribution : mu = %e,  sigma = %e \n',out_data.mu,out_data.sigma);
      fprintf(fid,' 75th percentile = %e \n', out_data.perc75);
      fprintf(fid,' median          = %e \n', out_data.median);
      fprintf(fid,' 25th percentile = %e \n', out_data.perc25);
      fprintf(fid,' skewness = %e, kurtosis = %e \n',out_data.skewness,out_data.kurtosis);
      fprintf(fid,'\n===================================================================\n');
    %
      fclose(fid);
    %
  end
  if(parFit.svd.svd==1)
      figure; semilogy((hist.Sc),'o'); hold on; semilogy((hist.Sr),'or'); hold off; legend('Sc','Sr');
      xlabel('Number of Singular Values');
      Titolo=['Threshold= ' num2str(parFit.svd.soglia,'%e')];
      title(Titolo); 
  end
  %
  catch1
  delete(h);
catch
  delete(h);
end
return;
