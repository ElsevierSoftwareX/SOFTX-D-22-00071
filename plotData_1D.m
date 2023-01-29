%#####################################################################################
%Name    : plotData_1D
%Purpose : To plot the measured data.
%Changes : 1.0 [28/07/2021] (vb)
%
%NOTES   :[28/07/2021] (vb) Derived from grafico_1D(x,T1,T2,method, FL_typeKernel)
%#####################################################################################
function plotData_1D(X,T, X_label, FileName, Titolo, Out_folder)
 %
 figure;
 plot(T,X,'LineWidth',1.5);
 grid on;
 set(gca,'FontSize',12,'fontweight','bold');
 axis([T(1) T(end) min(X) max(X)]);
 xlabel(X_label);
 ylabel('Intensity (a.u.)'); 
 title(Titolo);
 %
 OutputData=[T, X];
 %dlmwrite([Out_folder 'ArrayData.txt'],OutputData,'delimiter','\t','precision','%0.13e', 'newline', 'pc');
 dlmwrite([Out_folder FileName],OutputData,'delimiter','\t','precision',5, 'newline', 'pc');
 %
 return;
end
