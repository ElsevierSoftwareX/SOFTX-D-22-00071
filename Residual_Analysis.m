%###################################################################################################################################################
%NAME    : Residual_Analysis()
%PURPOSE : A posteriori descriptive statistics computation  
%
%VERSION : 1.0 [05/08/2021] (fz)
%          1.1 [08/08/2021] (vb) Now directly handles matrix X of computed data.
%          1.2 [18/08/2021] (vb) Changed the output plots.
%          1.3 [19/08/2021] (vb) Changed script name (Data_analysis => Residual_Analysis) and the list od formal parameters.
%          1.4 [12/10/2021] (vb) Changed the use of Flags and which statistic plots are showed (a few now are showed only under FL_debug)
%
%INPUT   : Res : Residual values
%                Res=kron(Kc,Kr)x-S  
%          X   : Matrix of computed data values
%                X=Kc*x*Kr'
%          S   : Matrix of data values
%OUTPUT  : out: structure   with fields:
%
%             out.sigma_R; % Residual Mean
%             out.mu_R;    % Residual variance
%             out.sigma;   % mean fitted Gaussian
%             out.mu;      % variance fitted Gaussian
%             out.perc75;  % 75th percentile
%             out.perc25;  % 25th percentile
%             out.skewness;
%             out.kurtosis;
%             out.pd;  output of fitdist matlab function
%
%NOTES   : A negative skewness value means the data is left skewed. 
%          The data has a larger peakedness than a normal distribution because the kurtosis value is greater than 3.
%
%#############################################################################################################################################################
%
%function out=Residual_Analysis(X, S, TreshPercent, Verbose)
function out=Residual_Analysis(X, S, parFL)
 %(vb-12/10/2021)
 TreshPercent= parFL.FL_Stat;
 Verbose=parFL.FL_Verbose;
 %
 if(Verbose)
   fprintf('---------------------------------------------------\n');
   fprintf('\n Compute descriptive statistics \n');
 end
 %
 %(vb 08/08/2021)
 % Measured and computed data discrepancy
 %-------------------------------------------------------------------------------------------------------------
 DataRes=(X-S);
 N_vals=numel(DataRes);
 RelDataResPerc=abs(DataRes./S)*100.0;
 size_i=size(RelDataResPerc,1);
 size_j=size(RelDataResPerc,2);
 %
 ImageDiscepancy=ones(size_i, size_j);
 for i =1:size_i
   for j=1: size_j
     if RelDataResPerc(i,j) > TreshPercent
        ImageDiscepancy(i,j)=0;
     end
   end
 end
 %(vb-12/10/2021)
 if parFL.FL_Debug
   %Relative Discrepancy map
   figure;%imagesc(RelDataResPerc);colormap gray;
   surf(RelDataResPerc);
   colorbar;
   xlabel('Matrix Data Discrepancy Columns'); 
   ylabel('Matrix Data Discrepancy Rows'); 
   title(['3D Relative Residual Map (%)']);
   %Binary relative discrepancy map
   figure;imagesc(ImageDiscepancy);colormap gray;
   xlabel('Matrix Data Discrepancy Columns'); 
   ylabel('Matrix Data Discrepancy Rows'); 
   title(['2D Binary Relative Residual Map, Black: Percentage Discrepancy > ', num2str(TreshPercent),' (%)']);
 end
 mu_R=sum(sum(DataRes))/N_vals;
 sigma_R=sum(sum((DataRes-mu_R).^2))/(N_vals-1);
 std_R=sqrt(sigma_R);
 figure;plot(1:N_vals,DataRes(:),'.',1:N_vals,mu_R*ones(N_vals,1),'-r');grid on;
 title(['1D Residual Plot, Mean Values  = ',num2str(mu_R),', std=',num2str(std_R)]);
  %--------------------------------------------------------------------------
 %Data statistics
  pd=fitdist(DataRes(:),'normal');
 sigma=pd.sigma;
 mu=pd.mu;
  figure; histfit(DataRes(:));grid on;title('Residual Histogram with computed Normal distribution fit');
 % compute the quantiles of the sample data.
 p = 0:0.25:1;
 y = quantile(DataRes(:),p);
 perc75=y(4);median_R=y(3);perc25=y(2);
 %(vb-12/10/2021)
 if parFL.FL_Debug
   %  Create a box plot to visualize the statistics.
   figure; boxplot(DataRes(:));title('Residual BoxPlot');grid on;
   %The box plot shows the 0.25, 0.5, and 0.75 quantiles. 
   %The long lower tail and plus signs show the lack of symmetry in the sample data values.
   %**************************************************************************
   % Compute the skewness and kurtosis of the data.
 end 
 sy = [skewness(DataRes(:)),kurtosis(DataRes(:))];
 %
 if(Verbose)
  fprintf('---------------------------------------------------\n');
  fprintf(' Fit Probability  Distribution    \n');
  fprintf(' parameter         confidence interval \n');
  fprintf('75th percentile : %e \n', perc75);
  fprintf('median: %e \n', median_R);
  fprintf('25th percentile : %e \n', perc25);
  fprintf('skewness = %e, kurtosis = %e \n',sy(1),sy(2));
  fprintf('---------------------------------------------------\n');
 end
 out.sigma_R=sigma_R;% Residual Mean computed with data measured and data reconstructed
 out.mu_R=mu_R;      % Residual sampling variance computed with data measured and data reconstructed
 out.sigma=sigma;    % mean fitted Gaussian
 out.mu=mu;          % var fitted Gaussian
 out.perc75=perc75;  % 75th percentile
 out.perc25=perc25;  % 25th percentile
 out.skewness=sy(1);
 out.kurtosis=sy(2);
 out.median=median_R;
 out.pd=pd;
 return;
end
