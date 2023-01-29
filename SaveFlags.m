%########################################################################################################################################
%NAME    : SaveFlags.m
%PURPOSE : Writes on file the values of the structured variabile parFL.
%VERSION : 1.0 [19/08/2021] (vb) Added FL_Stat. 
%          1.1 [03/02/2022] (vb) Added FL.FL_Lower_bound
%           
%NOTES   : 
%#########################################################################################################################################
function SaveFlags(data_folder,namefile,parFL,ExpName)
final_data=[data_folder namefile];fp=fopen(final_data,'w');
   fprintf(fp,'Flags for experiment %s \n%%\n',ExpName);
   fprintf(fp,'  FL_typeKernel         =%d\n',parFL.FL_typeKernel);
%    fprintf(fp,'FL_UPEN2D           =%d \n',parFL.FL_UPEN2D);
%    fprintf(fp,'FL_TIKHONOV         =%d \n',parFL.FL_TIKHONOV);
%    fprintf(fp,'FL_Stelar_Magriteck =%d \n',parFL.FL_Stelar_Magriteck);
   fprintf(fp,'  FL_InversionTimeLimits=%d\n',parFL.FL_InversionTimeLimits);
%    fprintf(fp,'FL_LoadMatrixB      =%d \n',parFL.FL_LoadMatrixB);
%    fprintf(fp,'FL_UseMatrixB       =%d \n',parFL.FL_UseMatrixB);
   fprintf(fp,'  FL_OutputData         =%d\n',parFL.FL_OutputData);
   fprintf(fp,'  FL_NoContour          =%d\n',parFL.FL_NoContour);
   fprintf(fp,'  FL_Verbose            =%d\n',parFL.FL_Verbose);
   fprintf(fp,'  FL_Debug              =%d\n',parFL.FL_Debug);
   fprintf(fp,'  FL_Amp_scale          =%d\n',parFL.FL_Amp_scale);
   fprintf(fp,'  FL_Lower_bound        =%d\n',parFL.FL_Lower_bound);
   fprintf(fp,'  FL_Scale_fact         =%d\n',parFL.FL_Scale_fact);
   fprintf(fp,'  FL_T1T2Filter         =%d\n',parFL.FL_T1T2Filter);
   fprintf(fp,'  FL_EraseCol           =%d\n',parFL.FL_EraseCol);
   fprintf(fp,'  FL_EraseRow           =%d\n',parFL.FL_EraseRow);
   fprintf(fp,'  FL_Offset             =%d\n',parFL.FL_Offset);
   fprintf(fp,'  FL_Stat               =%d\n',parFL.FL_Stat);
   fprintf(fp,'END');
   fclose(fp);
   return;
end