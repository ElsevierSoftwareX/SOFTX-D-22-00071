%#####################################################################################
%Name    : saveplot
%Changes : 1.1 [15/07/2021] (vb) Changed ...  .
%          1.2 [12/10/2021] (vb) Now doesn't save "*.fig".
%
%NOTES   :
%#####################################################################################
function saveplot(savefigpath,figname,GUIplot)
    plotname=[savefigpath '\' figname];
    fig1 = figure('visible','off');
    figAxes = axes(fig1);
    switch figname
                case '2D_Data_Input'
                    allChildren = GUIplot.XAxis.Parent.Children;
                    copyobj(allChildren, figAxes);
                    figAxes.XLabel = GUIplot.XLabel;
                    figAxes.YLabel = GUIplot.YLabel;
                    figAxes.XLim = GUIplot.XLim;
                    figAxes.YLim = GUIplot.YLim;
                    figAxes.XTick = GUIplot.XTick;
                    figAxes.XTickLabel = GUIplot.XTickLabel;
                    figAxes.YTick = GUIplot.YTick;
                    figAxes.YTickLabel = GUIplot.YTickLabel;
                    figAxes.DataAspectRatio = GUIplot.DataAspectRatio;
                    view(figAxes,55.5,21.6);
                    colorbar(figAxes);
                    grid(figAxes,'on'); 
                    box(figAxes,'on');                    
                case '2D_MAP' 
                    allChildren = GUIplot.XAxis.Parent.Children;
                    copyobj(allChildren, figAxes);
                    figAxes.XLabel = GUIplot.XLabel;
                    figAxes.YLabel = GUIplot.YLabel;
                    figAxes.XLim = GUIplot.XLim;
                    figAxes.YLim = GUIplot.YLim;
                    figAxes.XTick = GUIplot.XTick;
                    figAxes.XTickLabel = GUIplot.XTickLabel;
                    figAxes.YTick = GUIplot.YTick;
                    figAxes.YTickLabel = GUIplot.YTickLabel;
                    figAxes.DataAspectRatio = GUIplot.DataAspectRatio;
                    colorbar(figAxes);
                    grid(figAxes,'on'); box(figAxes,'on');
                case 'X_Projection'
                    allChildren = GUIplot.XAxis.Parent.Children;
                    copyobj(allChildren, figAxes);
                    figAxes.XLabel = GUIplot.XLabel;
                    figAxes.YLabel = GUIplot.YLabel;
                    figAxes.XScale = 'log';
                    grid(figAxes,'on'); box(figAxes,'on');
                case 'Y_Projection'
                    allChildren = GUIplot.XAxis.Parent.Children;
                    copyobj(allChildren, figAxes);
                    figAxes.XLabel = GUIplot.XLabel;
                    figAxes.YLabel = GUIplot.YLabel;
                    figAxes.XScale = 'log';
                    grid(figAxes,'on'); 
                    box(figAxes,'on');
                case 'Residuals_SVD'
                    allChildren = GUIplot.XAxis.Parent.Children;
                    copyobj(allChildren, figAxes);
                    figAxes.XLabel = GUIplot.XLabel;
                    figAxes.YLabel = GUIplot.YLabel;
                    figAxes.YScale = 'log';
                    legend('Sc','Sr');
    end
    %(vb-12/10/2021)
    %saveas(fig1,plotname);    % save it
    saveas(fig1,plotname,'png');
    close(fig1);             % clean up by closing it
    %delete(fig1);
return;
end