function SaveSetInputFile(data_folder,namefile,parFile,ExpName)
final_data=[data_folder namefile];
  fp=fopen(final_data,'w');
   fprintf(fp,'Filenames e parameters for experiment %s \n%%\n%% [File Data]\n',ExpName);
   fprintf(fp,'  filenamedata          =%s\n',parFile.filenamedata);
   fprintf(fp,'  filenameTimeX         =%s\n',parFile.filenameTimeX);
   fprintf(fp,'  filenameTimeY         =%s\n',parFile.filenameTimeY);
   %    fprintf(fp,'%% [Weight Matrix]\n  ');
   %    fprintf(fp,'FileNameMatrixB       =%s\n  ',parFile.FileNameMatrixB);
   fprintf(fp,'%% [Inversion Points]\n');   
   fprintf(fp,'  nx                    =%d\n',parFile.nx);
   fprintf(fp,'  ny                    =%d\n',parFile.ny);
   fprintf(fp,'%% [Inversion Time limits]\n');   
   fprintf(fp,'  T1min                 =%1.1e\n',parFile.T1min);
   fprintf(fp,'  T1max                 =%1.1e\n',parFile.T1max);
   fprintf(fp,'  T2min                 =%1.1e\n',parFile.T2min);
   fprintf(fp,'  T2max                 =%1.1e\n',parFile.T2max);
   fprintf(fp,'%%\nEND');
  fclose(fp);
  return;
end
