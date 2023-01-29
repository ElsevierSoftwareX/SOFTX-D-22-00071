% 
% MAIN for synthetic data creation.
%  Two peaks map described in 
%
%  Bortolotti V.; Landi G.; Zama F., 
%  2DNMR data inversion using locally adapted multi-penalty regularization, 
%  «COMPUTATIONAL GEOSCIENCES», 2021, 25, pp. 1215 - 1228
%
%  FZ 2021
close all
clear all
clc
FL_typeKernel=1;  
FL_UPEN2D=1;  
FL_TIKHONOV=0; 
DELTA=1.e-2;% Noise added
SOGLIA_T=10^(-2);%1.e-1;% 
B_mat=0;% No weights
par.B_mat=B_mat;
par.svd=0;%
delta=DELTA;
out_folder='./Synth_data_folder';mkdir(out_folder);
fprintf('------------------------------------------------------------------ \n')
fprintf(' Folder name %s \n',out_folder)
[Kc,Kr,s,gexact,par] =data2d_simul_2peaks_par(delta,2048, 128, 80,80);
    
Tau1=par.Tau1;Tau2=par.Tau2;scale_fact1=par.scaleTau1;
scale_fact2=par.scaleTau2;Amp_scale=par.scale_factAmpl; 
[nx,ny]=size(gexact);
% True relaxation map
Nome=[out_folder '/True_2Dmap.dat'];   
fprintf(' True map : %s \n',Nome)
dlmwrite(Nome,gexact,'delimiter','\t','precision','%0.13e')
%
% Relaxation Time channels in Y axis
%
Nome=[out_folder '/t_Y.dat'];
fprintf(' Relaxation Time channels in Y axis : %s \n',Nome)
fid=fopen(Nome,'w');
fprintf(fid,'%0.13e\n',par.Tau1);
fclose(fid);
%
% Relaxation Time channels in X axis
%
Nome=[out_folder '/t_X.dat'];
fprintf(' Relaxation Time channels in Y axis : %s  \n',Nome)
fid=fopen(Nome,'w');
fprintf(fid,'%0.13e\n',par.Tau2);
fclose(fid);


N=nx*ny;
normexact=norm((Kc*gexact*Kr'-s)/Amp_scale,'fro')^2;

[nx,ny]=size(gexact);

Nome=[out_folder '/s_ircpmg.dat'];   
fprintf(' 2D Relaxation data : %s  \n',Nome)

dlmwrite(Nome,s,'delimiter','\t','precision','%0.13e')     

fprintf('Scaled Residual Norm: %e \n',normexact)
fprintf('------------------------------------------------------------------ \n')
%%
grafico_1D(gexact,par.T1,par.T2,'True', 1)
grafico_2D(gexact,par.T1,par.T2,0, 1, 15, 'True')
grafico_3D(gexact,par.T1,par.T2, 1,  'True')

%==========================================================================
%  IR-CPMG synthetic data
%--------------------------------------------------------------------------
function [K1,K2,inputdata1,fmodel,par] = data2d_simul_2peaks_par(noise,BS,NBLK,nx,ny)
%
% IR-CPMG peak T1 > peak  T2
% 
%
 media_T1_1 = 2.5;% log10(Peak in T1) [-3 0] 
 sigma_T1_1 = 0.2;
 media_T2_1 = 2.;% log10(Peak in T2) [-3 0] 
 sigma_T2_1 = 0.2;
%
 media_T1_3 = 2.0;
 sigma_T1_3 = 0.2;
 media_T2_3 = 1.5;
 sigma_T2_3 = 0.5;
EDLY=250;  %us
PW90=3.8;  %us
ECHOSTEP=2*(EDLY+PW90)*10*1E-4; %it is in tenth of us.
%
Tau1 = logspace(-3,0.5,NBLK)'; 
Tau2max = 1;
Tau2 = linspace(0.001,Tau2max,BS)'; 
scale_factTau1=1E3; 
scale_factTau2=1E3; 
scale_factAmpl=1E4;
Tau1=scale_factTau1*Tau1; 
Tau2=scale_factTau2*Tau2;
Tau2(1)=ECHOSTEP;
for j=2:BS
  Tau2(j)=Tau2(j-1)+ECHOSTEP;
end
q1 = exp((1/(ny-1))*log(4*Tau1(end)/(0.25*Tau1(1))));
T1 = 0.25*Tau1(1)*q1.^(0:ny-1);T1_min=min(T1);T1_max=max(T1);
q2 = exp((1/(nx-1))*log(4*Tau2(end)/(0.25*Tau2(1))));
T2 = 0.25*Tau2(1)*q2.^(0:nx-1);T2_min=min(T2);T2_max=max(T2);
%------------------------------------------------------------------------------------
fprintf('Range T1 [%e, %e] log10 [%e,%e] \n',T1_min,T1_max,log10(T1_min),log10(T1_max));
fprintf('Range T2 [%e, %e] log10 [%e,%e] \n',T2_min,T2_max,log10(T2_min),log10(T2_max));
%
T1_exp = linspace(log10(T1_min),log10(T1_max),ny);%linspace(-3,0,ny);
T2_exp = linspace(log10(T2_min),log10(T2_max),nx);%linspace(-3,0,nx);
[fmodel1,PK1]=modello_componente(sigma_T1_1,T1_exp,media_T1_1,sigma_T2_1,T2_exp,media_T2_1);
[fmodel3,PK3]=modello_componente(sigma_T1_3,T1_exp,media_T1_3,sigma_T2_3,T2_exp,media_T2_3,pi/6);%Rotazione pi/6
fmodel = fmodel1+1.5*fmodel3 ; 
%--------------------------------------------------------------------------
fmodel = fmodel/max(max(fmodel))*128;
% cut-off
TT=1.E-10;
fmodel(fmodel<TT)=TT;
fmodel = fmodel-TT;

Kernel_1 = inline('1-2*exp(- Tau * (1./ TimeConst))','Tau','TimeConst');
K1 = Kernel_1 (Tau1,T1);
Kernel_2 = inline('exp( - Tau * (1./ TimeConst))','Tau','TimeConst');
K2 = Kernel_2(Tau2,T2);

% data
data = K1 * fmodel * K2';
dm = max(max(data));
data1 = data ./ dm;
fmodel = fmodel./dm; 
%
% Noise add
%
inputdata=addnoise(noise,data1);
inputdata1=inputdata.*scale_factAmpl;
fmodel=fmodel*scale_factAmpl;
inputdata=round(inputdata1);
%
DeltaQ = mean(abs(diff(T1)));
par.scale_factAmpl=scale_factAmpl;
par.scaleTau1=scale_factTau1;
par.scaleTau2=scale_factTau2;
par.DeltaQ=DeltaQ;
par.T1=T1;par.Tau1=Tau1;
par.T2=T2;par.Tau2=Tau2;
T1_t=par.T1;T2_t=par.T2;


end
%
%
function [y1]=addnoise(delta,y)
  rng('default');
  eta=randn(size(y));
  eta=eta/norm(eta(:));
  y1=y+sqrt(delta)*eta;
  fprintf(' \n Input Noise parameter delta=%e   \n Absolute difference between data and noisy data =%e \n',delta, norm(y1(:)-y(:))^2) 
end

function [fmodel,par]=modello_componente(sigma_T1,T1_exp,media_T1,sigma_T2,T2_exp,media_T2,theta)
    fmodelT1 = 1/(sigma_T1*((2*pi)^(1/2)))*exp(-((T1_exp-media_T1)./sigma_T1).^2/2);
    fac_T1=sum(fmodelT1);
    fmodelT2 = 1/(sigma_T2*((2*pi)^(1/2)))*exp(-((T2_exp-media_T2)./sigma_T2).^2/2);
    fac_T2=sum(fmodelT2);
 if nargin == 6
  theta=0;
 end
 Max_val=0;
  cy=media_T1;cx=media_T2;
  for i=1:length(T1_exp)
      Y=T1_exp(i)-cy;    
      for j=1:length(T2_exp)
          X=T2_exp(j)-cx;
          xr=cx+cos(theta)*X-sin(theta)*Y;
          yr=cy+sin(theta)*X+cos(theta)*Y;
          comp_T1=1/(sigma_T1*((2*pi)^(1/2)))*exp(-((xr-media_T1)./sigma_T1).^2/2);
          comp_T2=1/(sigma_T2*((2*pi)^(1/2)))*exp(-((yr-media_T2)./sigma_T2).^2/2);
          val=comp_T1*comp_T2/(fac_T2*fac_T1);fmodel(j,i)=val;
          if(val >= Max_val)
              Max_val=val;Peak=val;Y_pk=Y;X_pk=X;Row_pk=j;Col_pk=i;
          end
      end
  end
   par.Peak=Peak;par.Y=Y_pk;par.X=X_pk;par.I=Row_pk;par.J=Col_pk;
end





%==========================================================================
%  graphic functions
%--------------------------------------------------------------------------


function grafico_3D(X,tauv,tauh, FL_typeKernel,  Titolo)
    fig=figure;
    axes('FontSize',12);
    set(gcf,'Renderer','zbuffer');
    set(fig,'DoubleBuffer','on');
    set(gca,'NextPlot','replace','Visible','off')
    taulh = log10(tauh);
    sta = size(taulh);
    taulv = log10(tauv);
    stb = size(taulv);
    [XT1,XT2]=meshgrid(taulh,taulv);
    s=surf(XT1,XT2,X);alpha(s,'z');set(gca,'FontSize',12,'fontweight','bold')
    axis('tight')
    if (FL_typeKernel==1 || FL_typeKernel==2)
      xlabel('Log_{10}(T_2)  [T_2 in ms]'); 
      ylabel('Log_{10}(T_1)  [T_1 in ms]'); 
     elseif FL_typeKernel==3
      xlabel('Log_{10}(T_2)  [T_2 in ms]'); %
      ylabel('Log_{10} D (\mum^2/ms)'); %
     elseif FL_typeKernel==4
      xlabel('Log_{10}(T_{22})'); 
      ylabel('Log_{10}(T_{21})'); 
    end
    title(Titolo);
  return;
end   

function grafico_2D(X,tauv,tauh,Contur_or_Surf, FL_typeKernel, numline, Titolo)
    fig=figure;
    axes('FontSize',12);
    set(gcf,'Renderer','zbuffer');
    set(fig,'DoubleBuffer','on');
    set(gca,'NextPlot','replace','Visible','off')
    taulh = log10(tauh);
    sta = size(taulh);
    taulv = log10(tauv);
    stb = size(taulv);
    if Contur_or_Surf
      surf(taulh,taulv',X);hold on
      %plot(taulh,taulh,'--r');grid on;hold off
     else
      contour(taulh,taulv',X,numline,'Linewidth',1.3);hold on
      %plot(taulh,taulh,'--r','Linewidth',1.3);grid on;hold off
    end
    grid on
    caxis('auto');
    shading interp;
    axis([taulh(1),taulh(sta(2)),taulv(1),taulv(stb(2))]);
    colorbar;set(gca,'FontSize',12,'fontweight','bold')
    if (FL_typeKernel==1 || FL_typeKernel==2)
      xlabel('Log_{10}(T_2)  [T_2 in ms]'); 
      ylabel('Log_{10}(T_1)  [T_1 in ms]');
     elseif FL_typeKernel==3
      xlabel('Log_{10}(T_2)  [T_2 in ms]'); %
      ylabel('Log_{10} D (\mum^2/ms)'); %
     elseif FL_typeKernel==4
      xlabel('Log_{10}(T_{22})  [T_{22} in ms]'); 
      ylabel('Log_{10}(T_{21})  [T_{21} in ms]'); 
    end
    title(Titolo);
  return;
end   


function grafico_1D(x,T1,T2,metodo, FL_typeKernel)
 [nx,ny]=size(x);
 % peak
 [~,iy] = max(max(x));
 [~,ix] = max(max(x'));
 picco = x(ix,iy);
 M_picco=x(max(ix-5,1):min(ix+5,nx),max(iy-5,1):min(iy+5,ny)); 
 Perc=100*sum(M_picco(:))/sum(x(:));
 analisi_T1=sum(x,2);
 analisi_T2=sum(x,1);
 figure;
 semilogx(T1,analisi_T1,"LineWidth",2); 
 axis([T1(1) T1(end) min(analisi_T1) max(analisi_T1)]);
 set(gca,'FontSize',12,'fontweight','bold');grid on
 if (FL_typeKernel==1 || FL_typeKernel==2)
      xlabel('T_1  (ms)');
  elseif FL_typeKernel==3
      xlabel('D (\mum^2/ms)'); %
  elseif FL_typeKernel==4
      xlabel('T_{21}  (ms)'); 
 end
 ylabel('Intensity (a.u.)'); grid on;
 title(metodo);
 figure; semilogx(T2,analisi_T2,"LineWidth",2); 
 axis([T2(1) T2(end) min(analisi_T2) max(analisi_T2)]);
 set(gca,'FontSize',12,'fontweight','bold');grid on
 if (FL_typeKernel==1 || FL_typeKernel==2)
      xlabel('T_2  (ms)'); 
  elseif FL_typeKernel==3
      xlabel('T_2  (ms)'); %
  elseif FL_typeKernel==4
      xlabel('T_{22}  (ms)'); 
 end
 ylabel('Intensity (a.u.)');grid on;
 title(metodo);
 return;
end





                       
                       
