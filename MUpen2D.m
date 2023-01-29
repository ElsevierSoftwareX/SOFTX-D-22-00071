%#############################################################################################################################################################
%NAME    : MUpen2D.m
%PURPOSE : The main script for running the MUpen2D algorithm that uses both L1 and L2 regularizations
%VERSION : 1.0 [01/01/2020]
%          1.1 [11/07/2020] (vb) Now it works also with rectangular domain.
%                           Swapped the use of nx and ny between T1 and T2. T1 (or D)is on Y axis and T2 on X axis.
%                           Now nx is renamed ncol and ny is renamed nrow.
%                           In order to associate easily and clearly the variable T1, T2 (and D) to the XY axes of plots and to the matrix of data, 
%                           a "symbol protocol" must be established. Therefore:
%                           T1 corresponds to Y axis of the plot and to the rows of the matrix, therefore in the following it assumed that: (T1, Y, row);
%                           T2 corresponds to X axis of the plot and to the columns of the matrixt therefore in the following it assumed that: (T2, X, row);
%                                                    _____________________
%                                                   |                     | 
%                                                   |                     |
%                                                   |                     |
%                      Plot Y (matrix rows => T1, D)|                     |
%                                                   |                     |
%                                                   |                     |
%                                                   |                     |
%                                                   |                     |
%                                                   |_____________________|
%                                                  Plot X (matrix cols =>T2)
%          1.2 [20/07/2020] (vb) Added a priori knowledge constraint on T1-T2 map, "T2< 1.5 T1".
%                                Added new parameters read from filepar and the possibility to clean input data matrix.
%          1.3 [30/04/2021] (vb) Fixed a bug on size of matrix (Row and Col). Added a parameter(par.fista.weight) to weight 
%                                the influence of Upen and Fista.  
%          1.4 [05/05/2021] (vb) Fixed a bug on exit-loop controls. 
%          1.5 [11/05/2021] (fz)(gl) Fixed a bug of the Data_upen2d_fista()function. 
%          1.6 [19/05/2021] (vb) Now T1/T2 constraint work also on the function grad_proj_noreg().
%          1.7 [31/05/2021] (vb) Canceled the use of addpath(), not more supported by the R2020b compiler
%          1.8 [03/07/2021] (vb) Changed the structure of data par (in particular "FL_") and function to be more closed to the GUI approach. 
%                                Changed the name of the function and m-file.
%          1.9 [19/07/2021] (vb) Changed the management of parameter "par.fista.weight" and consmetic changes.
%          2.0 [21/07/2021] (vb) SVD singular values plot moved in the main.
%          2.0 [23/07/2021] (lb)(vb) Changed printout on file in function Data_upen2d_fista().
%          2.1 [26/07/2021] (vb) Now it is handled also the offset.
%          2.2 [28/07/2021] (vb) Changed the handling of compliance parameter beta_0.
%          2.3 [06/08/2021] (fz) Added statistic analysis and plots. 
%          2.4 [07/08/2021] (vb) Export as a parameter of hist structure the computed recontructed data matrix, "hist.X".
%          2.5 [19/08/2021] (vb) Now Upen loop uses only the tollerance for the distribution (TOL_x) to exit.
%                                Cosmetic changes to the code.
%          2.6 [20/08/2021] (vb) Changed the use of the T1T2filter, now it acts on the penatltie terms.
%          2.7 [25/08/2021] (vb) Implemented some change in functions Data_upen2d_fista() and Fista_Upen() to improve convergency. 
%                                In degub mode possibly Fista_Upen() divergency trend is monitored, and in case the routine is stopped.
%                                Now divergency is checked for 500 consecutive iterations.
%                                Fixed a bug in plotting grad_proj_noreg() matrix when Offset is on.
%                                Changed the list of formal parameters of Data_upen2d_fista().
%          3.0 [03/02/2022] (FZ) New parameter for map lower bound (allow for negative value). New kernel SR-SR.
%                                
%
%INPUT   : Test_Folder : folder containing data and parameters files.           
%OUTPUT  :            
%NOTES   : 
%#############################################################################################################################################################
function [x,T1,T2,hist]=MUpen2D(Test_folder, Data_folder, Out_folder, par, parFile, parFL, CommentTS, N_T1, N_T2, Tau1, Tau2, s)
 Amp_scale=parFL.FL_Amp_scale;
 scale_fact=parFL.FL_Scale_fact;
 %###################################### Set problem dimension ############################
 nRow=parFile.nx;
 nCol=parFile.ny;
 N=nCol*nRow;
 %
 %######################## Set times of the inversion channels ############################
 %Set times of the inversion channels. Two modalities: automatic setting or fixed setting.
 %Times in milliseconds
 %
 if(parFL.FL_InversionTimeLimits==1)
   Tau1=scale_fact*Tau1;
   Tau2=scale_fact*Tau2;
   q1 = exp((1/(nRow-1))*log(4*Tau1(end)/(0.25*Tau1(1))));
   T1 = 0.25*Tau1(1)*q1.^(0:nRow-1);
   q2 = exp((1/(nCol-1))*log(4*Tau2(end)/(0.25*Tau2(1))));
   T2 = 0.25*Tau2(1)*q2.^(0:nCol-1);
 else
   T1min=parFile.T1min;
   T1max=parFile.T1max;
   T2min=parFile.T2min;
   T2max=parFile.T2max;
   q1 = exp((1/(nRow-1))*log(T1max/T1min));
   T1 = T1min*q1.^(0:nRow-1);
   q2 = exp((1/(nCol-1))*log(T2max/T2min));
   T2 = T2min*q2.^(0:nCol-1);
 end
 %
 if parFL.FL_OutputData
  dlmwrite([Out_folder 'S.txt'],s,'delimiter','\t','precision',5, 'newline', 'pc');
  if(parFL.FL_typeKernel == 1||parFL.FL_typeKernel == 2)
    dlmwrite([Out_folder 't_1.txt'],Tau1,'delimiter','\t','precision',5, 'newline', 'pc');
    dlmwrite([Out_folder 't_2.txt'],Tau2,'delimiter','\t','precision',5, 'newline', 'pc');
    dlmwrite([Out_folder 'T1.txt'],T1','delimiter','\t','precision',5, 'newline', 'pc');
    dlmwrite([Out_folder 'T2.txt'],T2','delimiter','\t','precision',5, 'newline', 'pc');
  elseif(parFL.FL_typeKernel==4)
    dlmwrite([Out_folder 't_1.txt'],Tau1,'delimiter','\t','precision',5, 'newline', 'pc');
    dlmwrite([Out_folder 't_2.txt'],Tau2,'delimiter','\t','precision',5, 'newline', 'pc');
    dlmwrite([Out_folder 'T1.txt'],T1','delimiter','\t','precision',5, 'newline', 'pc');
    dlmwrite([Out_folder 'T2.txt'],T2','delimiter','\t','precision',5, 'newline', 'pc');
  end
   SaveAllPar(Out_folder, 'AllParUsed.txt',par,parFL, parFile,Test_folder); %(vb 21/07/2021)
   %SaveParFlags(Out_folder, 'ParFL.txt',parFL, Test_folder);
   %SaveParFile(Out_folder, 'ParFile.txt',parFile, Test_folder);
   %SaveParFit(Out_folder, 'ParFit.txt',par, Test_folder);
 end
 %
 %
 %############################# Set the Kernel #######################################
 if (parFL.FL_typeKernel==1|| parFL.FL_typeKernel==2) %SR or IR
    [Kernel_1,Kernel_2] = T1_T2_Kernel(parFL.FL_typeKernel);
  elseif(parFL.FL_typeKernel==3)  %D-T2
    [Kernel_1,Kernel_2] = D_T2_Kernel;
  elseif(parFL.FL_typeKernel==4)
    [Kernel_1,Kernel_2] = T2_T2_Kernel;
  elseif(parFL.FL_typeKernel==5)
    [Kernel_1,Kernel_2] = SR_SR_Kernel;%s=max(max(s))-s;%% DA TOGLIERE CON DATYI GIUSTI
 end
 Kc = Kernel_1 (Tau1,T1); 
 Kr = Kernel_2(Tau2,T2);
 %
 %add offset
 if(parFL.FL_Offset==1)
  Kc(:, end+1)=1; 
  Kr(:, end+1)=1;
 end
 % -------------------------------------------------------------------------
 % Inversion function
 % -------------------------------------------------------------------------
 %(vb-25/08/2021)
 %[x,LAMBDA,hist] = Data_upen2d_fista(Kc,Kr,s,par, T1, T2, parFL.FL_T1T2Filter, parFL.FL_typeKernel);
 [x,LAMBDA,hist] = Data_upen2d_fista(Kc,Kr,s,par, T1, T2, parFL);
 %
 % A priori knowledge constraints
 %
 %  if(parFL.FL_T1T2Filter&&(parFL.FL_typeKernel == 1||parFL.FL_typeKernel == 2))
 %    %A priori knowledge constraints if T2> 1.5 T1
 %    [size_T1]=size(T1,2);
 %    [size_T2]=size(T2,2);
 %    for i=1:size_T1
 %     for j=1:size_T2
 %      if T2(j)>2.0*T1(i)
 %        x(i,j)=0;  
 %      end
 %     end
 %     end
 %  end
 %
 Res_vec=Kc*x*Kr'-s;%
 %(vb-07/08/2021)
 hist.X=Kc*x*Kr';
 Res_final = norm(Res_vec,'fro')/norm(s,'fro');
 %
 %[vb - 11/10/2017] export computed distribution and parameters
 map_file=[Out_folder '2D_Distribution.txt'];
 fprintf('\n Final map File: %s \n',map_file);
 dlmwrite(map_file,x,'delimiter','\t','precision','%0.13e', 'newline', 'pc');
 % 
 final_data=[Out_folder 'Parameters.txt'];
 fprintf('\n Final Parameters file: %s \n\n',final_data);
 fp=fopen(final_data,'w');
   fprintf(fp,'--------------------------------------------------------------------------------------------------------- \n');
   %(vb-23/07/2021)
   %fprintf(fp,'MUpen2D Input Parameters \n upen_tol=%e,\n Projected Gradient Tol =%e \n Projected Newton Tol=%e \n Conjugate Gradient Tol =%e\n',...
   %        par.upen.tol,par.gpnr.tol,par.nwtp.tolrho, par.cgn2d.tol);
   fprintf(fp,'MUpen2D Input Parameters \n upen_tol=%e,\n Projected Gradient Tol =%e \n',par.upen.tol,par.gpnr.tol);    

   fprintf(fp,'SVD Threshold =%0.0e \n   Data size= %d x %d  \n',par.svd.soglia,hist.ssize(1),hist.ssize(2));
   fprintf(fp,'---------------------------------------------------------------------------------------------------------');
   fprintf(fp,'\r\n');
   fprintf(fp,'Number of Inversion channels:  horizontal %d, vertical  %d \n', nCol, nRow);   
   fprintf(fp,'Final Relative Residual Norm =%0.4e \r\n', Res_final);
   fprintf(fp,'Total MUpen2D Iterations = %d',numel(hist.it_int));
   fprintf(fp,'\r\n');
   fprintf(fp,'Total FISTA Iterations = %d ', sum(hist.it_int));
   fprintf(fp,'\r\n');
   fprintf(fp,'Computation Time = %4.5f s.',sum(hist.tempi));
   fprintf(fp,'\r\n');
   fprintf(fp,'---------------------------------------------------------------------------------------------------------\n');
   %fclose(fp);   
   % peak
   [~,ix] = max(max(x)); 
   [~,iy] = max(max(x')); 
   picco = x(iy,ix);
   %(vb)[07/07/2017]
   if (ix<=1) ix=2; end
   if (iy<=1) iy=2; end
   if (ix>=nCol) ix=nCol-1; end
   if (iy>=nRow) iy=nRow-1; end
   %M_picco=x(max(ix-5,1):min(ix+5,nCol),max(iy-5,1):min(iy+5,nRow)); 
   M_picco=x(max(iy-5,1):min(iy+5,nRow),max(ix-5,1):min(ix+5,nCol)); 
   Perc=100*sum(M_picco(:))/sum(x(:));
   if (parFL.FL_typeKernel==1||parFL.FL_typeKernel==2)
      fprintf(fp,'MUpen2D - T2=%0.2f T1=%0.2f peak=%0.2f  PercTot=%0.2f  \n',T2(ix),T1(iy),picco,Perc);
      fprintf(fp,'MUpen2D - (T2= %0.2f %0.2f %0.2f, T1 = %0.2f %0.2f %0.2f) \n',T2(ix-1),T2(ix),T2(ix+1),T1(iy-1),T1(iy),T1(iy+1));
    elseif (parFL.FL_typeKernel==4)
      fprintf(fp,'MUpen2D - T22=%0.2f T21=%0.2f peak=%0.2f  PercTot=%0.2f  \n',T2(ix),T1(iy),picco,Perc);
      fprintf(fp,'MUpen2D - (T22= %0.2f %0.2f %0.2f, T21 = %0.2f %0.2f %0.2f) \n',T2(ix-1),T2(ix),T2(ix+1),T1(iy-1),T1(iy),T1(iy+1));
   end
 fclose(fp);
 dlmwrite([Out_folder 'residual.txt'],Res_vec,'delimiter','\t','precision','%0.13e', 'newline', 'pc');
 % (vb 21/07/2021)
 %  grafico_1D(x,T1,T2, '1D Distribution', parFL.FL_typeKernel, Out_folder);
 %  grafico_2D(x,T1,T2,0 ,parFL.FL_typeKernel, '2D Map');
 %  grafico_2D(x,T1,T2,parFL.FL_NoContour ,parFL.FL_typeKernel, '2D Map');
 %  grafico_3D(x,T1,T2,parFL.FL_typeKernel, '3D Distribution');
 %
 return;
end
%
%
%########################################################################################################################################
%NAME    : Data_upen2d_fista()
%PURPOSE : Implementation of the multi-penalty  inversion algorithm proposed in
%           Bortolotti V.; Landi G.; Zama F., 2DNMR data inversion using locally adapted multi-penalty regularization,
%           «COMPUTATIONAL GEOSCIENCES», 2021, 25, pp. 1215 - 1228
%
%VERSION : 1.0
%          1.1 [25/08/2021](vb)Changed the list of formal parameters, now it used the parFL structure.
%                              Parameter lambda is computed in accordance with eq. 21 (alfa), so that the residuals from
%                              regularization L1 and L2 are comparable.
%                              Denormalization of parameter beta_0 performed in the loop.
%
%INPUT   : Kc: kernel T1
%          Kr: kernel T2
%          s : data
%          par : structure for algorithm parameters (see manual)
%          T1, T2 : 
%          parFL: structure with flag parameters (see manual)
%OUTPUT  : x,
%          ck:
%          hist : structure containing significant computed values
%                fields:
%                Sc: singular values of Kc
%                Sr: singular values of Kr
%                res_vec: final residual 
%                ssize: size of data s
%                res: residual norm at outer iteration 
%                res_int: cell array, residual norm inside FISTA iterations
%                obj: cell array, value of the objective function inside FISTA iterations
%                tempi: time vector of  FISTA steps
%                it_int: vector: number of FISTA iterations
%                lambda_inner: vector values of L1 regularization parameter
%                it_cg: number of Initial Projected Gradient iterations
%                ck: final local L2 regularization parameters
%
%NOTES   :
%
%########################################################################################################################################
%function [x,ck,hist]=Data_upen2d_fista(Kc,Kr,s,par, T1, T2, FL_T1T2Filter, FL_typeKernel)
function [x,ck,hist]=Data_upen2d_fista(Kc,Kr,s,par, T1, T2, parFL)

  scale_fact = par.scale_fact;
  TOL_x      = par.upen.tol;
  TOL_res    = par.upen.tol_res; %not used
  Kmax       = par.upen.iter;
  if isfield(par,'Amp_scale')
      Amp_scale=par.Amp_scale;
    else
      Amp_scale=1;
  end
  %(FZ)[05/08/2021]  L2-L1 Balancing weights
  %      if    0 <= par.fista.weight  <= 1
  %                       L1 weight = par.fista.weight
  %                       L2 weight = 1- par.fista.weight
  %      else
  %           L1 weight = L2 weight =1  (no weighting factor)
  %      end
  %
  if (par.fista.weight>=0)&&(par.fista.weight<=1)
    WeightFisUpenL1=par.fista.weight;
    WeightFisUpenL2=1.0-par.fista.weight;
   else
    WeightFisUpenL1=1.0;
    WeightFisUpenL2=1.0;
  end
  %--------------------------------------------------------------------------
  % 
  [Uc,Sc,Vc]=svd(Kc); 
  [Ur,Sr,Vr]=svd(Kr); 
  %
  %
  L=Sc(1,1)*Sr(1,1);
  L=L*L;
  % TSVD filer/projection
  if par.svd.svd
    %--------------------------------------------------------------------------
    soglia=par.svd.soglia;
    Sc=diag(Sc);
    hist.Sc=Sc;
    Sr=diag(Sr);
    hist.Sr=Sr;
    if soglia < min(Sc)
     nc=length(Sc);
    else
     nc=find(Sc<=soglia,1);
    end
    if soglia < min(Sr)
      nr=length(Sr);
    else
      nr=find(Sr<=soglia,1);
    end
    %--------------------------------------------------------------------------
    Uc=Uc(:,1:nc); Vc=Vc(:,1:nc); Sc=Sc(1:nc);
    Ur=Ur(:,1:nr); Vr=Vr(:,1:nr); Sr=Sr(1:nr);
    s=Uc'*s*Ur;
    Kc=Uc'*Kc;
    Kr=Ur'*Kr;
 end
 %--------------------------------------------------------------------------
 [N_T1,N_T2]=size(s); 
 NN = (N_T1*N_T2)+1;% 
 nRow = size(Kc,2); 
 nCol = size(Kr,2);
 ck=ones(nRow, nCol)*1.0E10;
 [L1nx,L1ny,L2] = get_diff(nRow,nCol);
 % Set up parameters
 beta_00 = par.upen.beta00;
 beta_p  = beta_00*par.upen.beta_p; 
 beta_c  = beta_00*par.upen.beta_c;
 %(vb 28/07/2021)
 beta_0  =par.upen.beta0;
  %beta_0  = beta_00*par.upen.beta0;
 %
 [x0, iter]=grad_proj_noreg(Kc,Kr,s,parFL.FL_Lower_bound, zeros(nRow,nCol), par);
 
  %(fz-03-02-2022) do not apply the filter to general data (not NMR)
  if parFL.FL_Lower_bound==0
  %(vb-18/07/2020) application of a filter to reduce border effects.
   x0(:, 1:3)=0;
   x0(:, (end-3):end)=0;
   x0(1:3, :)=0;
   x0((end-3):end,:)=0;
  end
 %[19/05/2021]
 if(parFL.FL_T1T2Filter&&(parFL.FL_typeKernel == 1||parFL.FL_typeKernel == 2))
   %A priori knowledge constraints if T2> 1.5 T1
   [size_T1]=size(T1,2);
   [size_T2]=size(T2,2);
   for l=1:size_T1
    for j=1:size_T2
     if T2(j)>2.0*T1(l)
       x(l,j)=0;  
     end
    end
   end
 end
 %[vb -25/08/2021]
 if parFL.FL_Debug
   x1=x0;
   if(parFL.FL_Offset==1)
     x1(:, end)=[];
     x1(end,:)=[];
   end
   grafico_2D(x1,T1,T2,30 ,par.Kernel, '2D Map Grad Projec');
 end
 %
 x = x0;
 Rsqrd = norm((Kc*x*Kr'-s),'fro')^2; 
 % 
 c = reshape(L2*(x(:)),nRow,nCol); 
 c = ordfilt2(abs(c),9,ones(3));
 % 
 px = L1nx*(x(:)); py = L1ny*(x(:)); 
 v = sqrt(px.^2+py.^2);
 p = reshape(v,nRow,nCol);
 p = ordfilt2(p,9,ones(3));
 %(vb-28/07/2021)
 beta_0 = par.upen.beta0*Rsqrd/NN; % it must be a threshold const parameter in the computation of ck
 ck = Rsqrd./(NN*(beta_0+beta_p*p.^2+beta_c*c.^2)); 
 %(vb)[18/07/2021]
 %ck=(1-WeightFisUpen)*ck;
 ck=WeightFisUpenL2*ck;
 %(vb-25/08/2021)
 %lambda=WeightFisUpen*(Rsqrd/(norm(x(:),1)*NN));
 lambda=WeightFisUpenL1*(Rsqrd/(norm(x(:),1)));
 %
 continua = 1; 
 i=1;
 rest_diff=0;
 %cond_res=1;
 %
 %X_init=zeros(size(x0));
 %X_init=x0;
 %x=X_init;
 %x =zeros(size(x0));
 x=x0;
 %
 while continua
    xold = x; 
    [x, par_out] = Fista_Upen(x,s,Kr,Kc,ck,L2,par,L,lambda,parFL.FL_Lower_bound);
    res_int{i}=par_out.nres;
    tempi(i)=par_out.times;
    fobj{i}=par_out.objective;
    it_inter(i)=numel(par_out.nres);
    x = reshape(x,nRow,nCol);
    c = reshape(L2*(x(:)),nRow,nCol); 
    c = ordfilt2(abs(c),9,ones(3));
    % 
    px = L1nx*(x(:)); py = L1ny*(x(:)); 
    v = sqrt(px.^2+py.^2);
    p = reshape(v,nRow,nCol);
    p = ordfilt2(p,9,ones(3));
    % output
    res(i) = norm((Kc*x*Kr'-s)/Amp_scale,'fro');
    %(vb-24/08/2021)
    %distribution_diff=norm(x-xold,'fro')/norm(x,'fro');
    distribution_diff=norm(x-xold,'fro')/norm(xold,'fro');
    cond_err=distribution_diff>TOL_x;
    if (par.Verbose)
      fprintf('Data_upen2d_fista iter = %d, distribution differenze=%3.3g\n',i,distribution_diff);
    end
     continua = cond_err && i<Kmax;
    if continua
      Rsqrd = norm((Kc*x*Kr'-s),'fro')^2; 
      beta_0 =par.upen.beta0*Rsqrd/NN; 
      ck =Rsqrd./(NN*(beta_0+beta_p*p.^2+beta_c*c.^2)); 
      ck=WeightFisUpenL2*ck;
      lambda=WeightFisUpenL1*(Rsqrd/(norm(x(:),1)));
      %
      lambda_inner(i)=lambda;
      i = i+1;
    elseif(i==1)  %(vb)[05/05/2021]
      lambda_inner(1)=lambda;   %if the loop ends at first step (i==1) then lambda_inner is undefined.
    end
    %(vb 20/08/2021)
    if(parFL.FL_T1T2Filter&&(parFL.FL_typeKernel == 1||parFL.FL_typeKernel == 2))
      for l=1:size_T1
        for j=1:size_T2
         if T2(j)>2.0*T1(l)
           ck(l,j)=1.0E10;  
         end
        end
     end
    end
 end
 %---------------------
 % FZ added 03/08/2021
 res_vec=Kc*x*Kr'-s;
 hist.res_vec=res_vec;
 %----------------------
 hist.ssize=[N_T1,N_T2];
 hist.res = res;
 hist.res_int=res_int;
 hist.obj= fobj;
 hist.tempi=tempi;
 hist.it_int=it_inter;
 hist.lambda_inner=lambda_inner;
 hist.it_cg=iter;
 hist.ck=ck;
 return;
end
%
%
%########################################################################################################################################
%NAME    : Fista_Upen()
%PURPOSE : Implementation of the adapted FISTA algorithm as reported in
%          Algorithm FISTA_STEP reported in
%           Bortolotti V.; Landi G.; Zama F., 2DNMR data inversion using locally adapted multi-penalty regularization,
%           «COMPUTATIONAL GEOSCIENCES», 2021, 25, pp. 1215 - 1228
%
%VERSION : 1.01 24/08/2021 (vb) Changed wheight for terms temp_Upen and temp_Fista (1.5 and 0.5, respectively), the total fit error must be equal to 2 epsilon^2 (Lemma 2.1).
%                               Added additional printout in verbose mode.
%
%INPUT   : x,b,Kr,Kc,ck,L2,par,L,tau    
%
%OUTPUT  : x: solution f^(k+1) in  Algorithm FISTA_STEP.
%          par_out: structure of output computed values. Fields:
%          par_out.it_int: number of FISTA iterations 
%                         (index j in Algorithm FISTA_STEP)
%          par_out.objective : values of the objective function at each iteration;
%          par_out.times : computation time
%          par_out.nres : residual norm at each iteration
%NOTES   : 
%
%#############################################################################################################################################################
function [x, par_out] = Fista_Upen(x,b,Kr,Kc,ck,L2,par,L,tau,lb)
 par_out=struct;
 tolerance=par.fista.tol;
 maxiters=par.fista.maxiter;
 Amp_scale=par.Amp_scale;
 %
 [~,nCol]=size(Kr); 
 [~,nRow]=size(Kc);
 b=[b(:); zeros(nCol*nRow,1)];   
 T_reg=sqrt(ck);
 x=x(:);
 y = x;
 t = 1;
 L_A=L;
 Lc=64*max(max(abs(ck)));
 L=L_A+Lc;
 %
 tic
 temp_Upen=1.5*norm(A(x,Kr,Kc,L2,T_reg)-b)^2;
 %(vb-24/08/2021)
 %temp_Fista=tau*sum(abs(x));
 temp_Fista=0.5*tau*sum(abs(x));
 objective(1)=temp_Upen+temp_Fista;
 %
 if(par.Verbose)
    fprintf('\n\n   Fista_Upen iter = %d, L =%3.3g, Rel_obj_diff= 0, obj_Upen = %10.8g, obj_Fista=%10.8g\n',1,L,temp_Upen,temp_Fista);
 end
 %criterion of stop if FISTA relative residual increase for more then 100 consecutive times.
 Conv=0;
 continua=1;
 criterion_old=0;
 objective(1)=1E10;
 k=1;
 %j=0;
 while continua
    k=k+1;
    x_old = x;
    t_old = t;
    y1=A(y,Kr,Kc,L2,T_reg) - b; 
    y = y - (1/L)*AT(y1,Kr,Kc,L2,T_reg );
    x = soft(y,tau/L);
    t = 0.5*(1 + sqrt(1+4*t_old^2));
    y = x + ((t_old - 1)/t )*(x - x_old);
    nres(k)=norm(y1/Amp_scale);
    %
    temp_Upen=1.5*norm(A(x,Kr,Kc,L2,T_reg)-b)^2;
    temp_Fista=0.5*tau*sum(abs(x));
    objective(k) =temp_Upen+temp_Fista;
    %
    criterion = abs(objective(k)-objective(k-1))/objective(k-1);
    continua= criterion > tolerance && k <= maxiters;
    if(criterion_old<criterion&&par.Debug)% &&k>1)
       Conv=Conv+1;
       if Conv==500 
           %continua=0;
           if (par.Verbose)
              fprintf('\n   Fista_Upen iter = %d, Exit with Conv >%d \n',k, Conv);
           end
           break;
       end
    else
       Conv=0;
    end
    if (par.Verbose)
       fprintf('   Fista_Upen iter = %d, Rel_obj_diff=%5.5g, obj_Upen = %14.10g, obj_Fista=%14.10g \n',k, criterion, temp_Upen,temp_Fista);
    end
    %
    criterion_old=criterion;
 %
 end
 times=toc;
 par_out.it_int=k;
 par_out.objective= objective;
 par_out.times= times;
 par_out.nres=nres;
 %
 return;
end
%
%
%########################################################################################################################################
% NAME    : soft
% PURPOSE : Implementation of soft thresholding
%
%
%INPUT   : x vector  z^(j) in [*] p. 1218
%        : T threshold value
%OUTPUT  : y vector  f^(j)  in  [*] p. 1218      
%NOTES   : 
%         [*]           Bortolotti V.; Landi G.; Zama F., 2DNMR data inversion using locally adapted multi-penalty regularization,
%                      «COMPUTATIONAL GEOSCIENCES», 2021, 25, pp. 1215 - 1228
%#############################################################################################################################################################
function y = soft(x,T)
  y=sign(x-T).*(max((abs(x)-T),0));
  return;
end
%
%
%########################################################################################################################################
% NAME    : A
% PURPOSE : Implementation of matrix-vector product 
%             -                     -
%      y =   | kron(Kc,Kr),|   x 
%            | L*theta2    |
%             -           
%
%INPUT   :x,Kr,Kc,L2,theta2
%OUTPUT  : y     
%NOTES   : 
%#############################################################################################################################################################
function y=A(x,Kr,Kc,L2,theta2)
 [~,nRow]=size(Kc); 
 [~,nCol]=size(Kr);
 x_reshaped=reshape(x,nRow,nCol);
 y = Kc*x_reshaped*Kr'; 
 y=y(:);
 y=[y;theta2(:).*(L2*x)];
 return;
end
%
%########################################################################################################################################
% NAME    : AT()
% PURPOSE : Implementation of matrix transpose-vector product 
% 
%             -                          -
%      y =   | kron(Kc,Kr)^T, 0           |   x 
%            | 0,           (L*theta2)^T  |
%             -                          -
%
%
%INPUT   :x,Kr,Kc,L2,theta2
%OUTPUT  : y     
%NOTES   : 
%#############################################################################################################################################################
function y=AT(x,Kr,Kc,L2,theta2)
 % x  un vettore
 %  [N1,nx]=size(Kr);
 %  [N2,ny]=size(Kc);
 %  VV=reshape(x(1:N1*N2),N2,N1);
 %  y=Kc'*VV*Kr; y=y(:);  %% Correzione FZ 10/10/17
 %  y=y+theta2(:).*(L2*x(N1*N2+1:nx*ny+N1*N2));
 [N2,nRow]=size(Kc);
 [N1,nCol]=size(Kr);
 VV=reshape(x(1:N1*N2),N2,N1);
 y=Kc'*VV*Kr; y=y(:);  %% Correzione FZ 10/10/17
 y=y+theta2(:).*(L2*x(N1*N2+1:nRow*nCol+N1*N2));
 return;
end
%
%
%############################################################################################################################################################
% NAME    : grad_proj_noreg()
% PURPOSE : Implementation of projected gradient algorithm to solve
%
%   min_x>0 || A x - s ||,   A=kron(Kc,Kr)
%
%
%INPUT   : Kc,Kr,s, lb, x0, par
%OUTPUT  : x, k, norma_grad     
%NOTES   : 
%#############################################################################################################################################################
function [x, k, norma_grad]=grad_proj_noreg(Kc,Kr,s, lb, x0, par)
 %
 % set up tolerances
 tol=par.gpnr.tol;
 maxk=par.gpnr.maxiter;
 %STEP 1 (starting guess)
 x = max(x0,lb);
 %[nx,ny] = size(x0);
 [nRow,nCol] = size(x0);
 alpha_min=1.0E-10; 
 alpha_max=1.0E10; 
 alpha=1;        
 % Gradient of the objective function
 temp = Kc*x*Kr'-s; 
 res = temp;
 grad=Kc'*temp*Kr;   %A'(Ax-b) 
 norma_res(1)=norm(res(:)/par.Amp_scale); 
 k=1; 
 continua = 1;
 while continua
    %STEP 2 (Projection)
    d=max(x-alpha*grad,lb)-x;
    %STEP 3
    temp = Kc*d*Kr'; 
    temp=Kc'*temp*Kr;
    Ad = temp;
    if norm(Ad(:))>eps*norm(d(:))
        lambda=min(-(grad(:)'*d(:))/(d(:)'*Ad(:)), 1);     
      else
        lambda=1;
    end
    x=x+lambda*d;
    grad=grad+lambda*Ad;
    res = Kc*x*Kr'-s;
    %STEP 4 
    if norm(Ad(:))>eps*norm(d(:))
        if mod(k, 6)<3    
            alpha=(d(:)'*Ad(:))/(Ad(:)'*Ad(:));
          else
            alpha=(d(:)'*d(:))/(d(:)'*Ad(:));
        end
        alpha=max(alpha_min,min(alpha_max, alpha));
      else
        alpha=alpha_max;
    end
    k=k+1;
    norma_grad(k)=norm(grad(:));  
    norma_res(k)=norm(res(:)/par.Amp_scale); 
    Diff_res= abs(norma_res(k)-norma_res(k-1));
    continua = k<maxk && Diff_res>=tol;
 end 
 if (k >= maxk && par.Verbose)
   fprintf('*** Exit max iter grad proj noreg Diff res = %e  k=%d \n',Diff_res, k)
 end
 return;
end
