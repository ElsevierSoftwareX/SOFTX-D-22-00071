%########################################################################################################################################
%NAME    : SaveParFit.m
%PURPOSE : Writes on file the values of the structured variabile parFit.
%VERSION : 1.0 [21/07/2020] (vb) cosmetic changes.
%
%NOTES   : 
%#########################################################################################################################################
function SaveParFit(data_folder,namefile,par,ExpName)
 final_data=[data_folder namefile];
 ExpName=['"' ExpName '"'];
 fp=fopen(final_data,'w');
   fprintf(fp,'Fitting parameters for the experiment in folder %s \n\n',ExpName);
   fprintf(fp,'%%[Projected Gradient]\n');
   fprintf(fp,'  par.gpnr.tol          =%1.1e\n',par.gpnr.tol);
   fprintf(fp,'  par.gpnr.maxiter      =%1.1e\n',par.gpnr.maxiter);
   fprintf(fp,'%%[Projected Newton]\n  ');   
   fprintf(fp,'  par.nwtp.maxiter      =%d\n',par.nwtp.maxiter);
   fprintf(fp,'  par.nwtp.tolrho       =%1.1e\n',par.nwtp.tolrho);
   %    fprintf(fp,'par.nwtp.psi          =%1.1e\n  ',par.nwtp.psi);
   %    fprintf(fp,'par.nwtp.maxarm       =%d\n',par.nwtp.maxarm);
   %    fprintf(fp,'par.nwtp.eta          =%1.1e\n',par.nwtp.eta);
   fprintf(fp,'%%[Conjugate Gradient]\n  ');   
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
   fprintf(fp,'%%\nEND');
 fclose(fp);
 return;
end
