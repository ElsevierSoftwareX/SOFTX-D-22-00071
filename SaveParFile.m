%########################################################################################################################################
%NAME    : SaveParFile.m
%PURPOSE : Writes on file the values of the structured variabile parFile.
%VERSION : 1.0 [21/07/2020] (vb) cosmetic changes.
%NOTES   : 
%#########################################################################################################################################
function SaveParFile(data_folder,namefile,parFile,ExpName)
  final_data=[data_folder namefile];
  ExpName=['"' ExpName '"'];
  fp=fopen(final_data,'w');
    fprintf(fp,'Filenames and parameters for the experiment in folder %s \n\n',ExpName);
    fprintf(fp,'%% \n');
    fprintf(fp,'%% [File Data] \n');
    fprintf(fp,'  filenamedata          =%s\n',parFile.filenamedata);
    fprintf(fp,'  filenameTimeX         =%s\n',parFile.filenameTimeX);
    fprintf(fp,'  filenameTimeY         =%s\n',parFile.filenameTimeY);
    fprintf(fp,'%%\n');
    fprintf(fp,'%% [Inversion Points] \n');
    fprintf(fp,'  nx                    =%d\n',parFile.nx);
    fprintf(fp,'  ny                    =%d\n',parFile.ny);
    fprintf(fp,'%% [Inversion Time limits] \n');
    fprintf(fp,'  T1min                 =%1.1e\n',parFile.T1min);
    fprintf(fp,'  T1max                 =%1.1e\n',parFile.T1max);
    fprintf(fp,'  T2min                 =%1.1e\n',parFile.T2min);
    fprintf(fp,'  T2max                 =%1.1e\n',parFile.T2max);
    fprintf(fp,'%% \n');
    fprintf(fp,'END');
  fclose(fp);
  return;
end