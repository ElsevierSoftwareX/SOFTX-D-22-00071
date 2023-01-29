%########################################################################################################################################
%NAME    : SaveParFit.m
%PURPOSE : Writes on file the values of the structured variabile parFit.
%VERSION : 1.0 [21/07/2020] (vb) cosmetic changes.
%          1.1 [19/08/2021] (vb) Added FL_Stat. 
%          1.2 [03/02/2022] (vb) Added FL.FL_Lower_bound
%NOTES   : 
%#########################################################################################################################################
function SaveAllPar(data_folder,namefile,par, parFL, parFile,ExpName)
 final_data=[data_folder namefile];
 ExpName=['"' ExpName '"'];
 fp=fopen(final_data,'w');
   fprintf(fp,'All parameters used for the experiment in folder %s \n\n',ExpName);
   %
   fprintf(fp,'Used Fitting parameters\n');
   fprintf(fp,'%%[Projected Gradient]\n');
   fprintf(fp,'  par.gpnr.tol          =%1.1e\n',par.gpnr.tol);
   fprintf(fp,'  par.gpnr.maxiter      =%1.1e\n',par.gpnr.maxiter);
   fprintf(fp,'%%[Projected Newton]\n');   
   fprintf(fp,'  par.nwtp.maxiter      =%d\n',par.nwtp.maxiter);
   fprintf(fp,'  par.nwtp.tolrho       =%1.1e\n',par.nwtp.tolrho);
   %    fprintf(fp,'par.nwtp.psi          =%1.1e\n  ',par.nwtp.psi);
   %    fprintf(fp,'par.nwtp.maxarm       =%d\n',par.nwtp.maxarm);
   %    fprintf(fp,'par.nwtp.eta          =%1.1e\n',par.nwtp.eta);
   fprintf(fp,'%%[Conjugate Gradient]\n');   
   fprintf(fp,'  par.cgn2d.tol         =%1.1e\n',par.cgn2d.tol);
   fprintf(fp,'  par.cgn2d.maxiter     =%1.1e\n',par.cgn2d.maxiter);
   fprintf(fp,'%%\n%%[VARIOUS]\n%%[SVD]\n');   
   fprintf(fp,'  par.svd.svd           =%d\n',par.svd.svd);
   fprintf(fp,'  par.svd.soglia        =%1.1e\n',par.svd.soglia);
   fprintf(fp,'%%\n%%[UPEN]\n');   
   fprintf(fp,'  par.upen.tol          =%1.1e\n',par.upen.tol);
   fprintf(fp,'  par.upen.iter         =%d\n',par.upen.iter);
   fprintf(fp,'  par.upen.tol_res      =%1.1e\n',par.upen.tol_res);
   fprintf(fp,'  par.upen.beta00       =%1.1e\n',par.upen.beta00);
   fprintf(fp,'  par.upen.beta0        =%1.1e\n',par.upen.beta0);
   fprintf(fp,'  par.upen.beta_p       =%1.1e\n',par.upen.beta_p); 
   fprintf(fp,'  par.upen.beta_c       =%1.1e\n',par.upen.beta_c);
   fprintf(fp,'%%\n%%[FISTA]\n'); 
   fprintf(fp,'  par.fista.maxiter     =%1.1e\n',par.fista.maxiter);
   fprintf(fp,'  par.fista.tol         =%1.1e\n',par.fista.tol); 
   fprintf(fp,'  par.fista.weight      =%1.1e\n',par.fista.weight);
   fprintf(fp,'%%\n%%[TIKHONOV]\n');
   fprintf(fp,'  par.tikh.lambda       =%1.1e\n',par.tikh.lambda);   
   fprintf(fp,'%%\nEND Fit\n\n\n');
   %
   fprintf(fp,'Used Flags parameters\n');
   fprintf(fp,'  FL_typeKernel         =%d\n',parFL.FL_typeKernel);
   fprintf(fp,'  FL_InversionTimeLimits=%d\n',parFL.FL_InversionTimeLimits);
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
   fprintf(fp,'END Flags\n\n\n');
   %
   fprintf(fp,'Used Filenames and parameters\n');
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
   fprintf(fp,'%%\n');
   fprintf(fp,'END Filenames and parameters');  
 fclose(fp);
 return;
end
