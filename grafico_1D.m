%#####################################################################################
%Name    : grafico_1D
%Changes : 1.1 [03/05/2019] (vb) Changed the name of the script.
%          1.2 [25/07/2019] (fz) increased font  and line size.
%          1.3 [18/07/2020] (vb) Swapped the use of nx and ny between T1 and T2. 
%                                T1 (or D)is on Y axis and T2 on X axis of the chart.
%          1.4 [20/05/2021] (vb) Now also exports on file the T1 and T2 distributions. 
%
%NOTES   :[13/12/2017] (vb) Derived from grafico_1(x,T1,T2,method, FL_typeKernel)
%#####################################################################################
function grafico_1D(x,T1,T2,metodo, FL_typeKernel, Out_folder)
 [ny,nx]=size(x);
 % peak
 [~,ix] = max(max(x));
 [~,iy] = max(max(x'));
 picco = x(iy,ix);
 M_picco=x(max(iy-5,1):min(iy+5,ny),max(ix-5,1):min(ix+5,nx)); 
 Perc=100*sum(M_picco(:))/sum(x(:));
 analisi_T1=sum(x,2);
 analisi_T2=sum(x,1);
 figure;
 semilogx(T1,analisi_T1,'LineWidth',1.5);
 grid on;
 set(gca,'FontSize',12,'fontweight','bold');
 axis([T1(1) T1(end) min(analisi_T1) max(analisi_T1)]);
 if (FL_typeKernel==1 || FL_typeKernel==2)
      xlabel('T_1  (ms)');
  elseif FL_typeKernel==3
      xlabel('D (\mum^2/ms)'); %
  elseif FL_typeKernel==4
      xlabel('T_{21}  (ms)'); 
 end
 ylabel('Intensity (a.u.)'); 
 title(metodo);
 figure; semilogx(T2,analisi_T2,'LineWidth',1.5); 
 set(gca,'FontSize',12,'fontweight','bold')
 axis([T2(1) T2(end) min(analisi_T2) max(analisi_T2)]);grid on
  if (FL_typeKernel==1 || FL_typeKernel==2)
      xlabel('T_2  (ms)'); 
  elseif FL_typeKernel==3
      xlabel('T_2  (ms)'); %
  elseif FL_typeKernel==4
      xlabel('T_{22}  (ms)'); 
 end
 ylabel('Intensity (a.u.)');
 title(metodo);
 %
 dlmwrite([Out_folder 'T1_distribution.txt'],analisi_T1,'delimiter','\t','precision','%0.13e', 'newline', 'pc');
 dlmwrite([Out_folder 'T2_distribution.txt'],analisi_T2,'delimiter','\t','precision','%0.13e', 'newline', 'pc');
 %
  return;
end
