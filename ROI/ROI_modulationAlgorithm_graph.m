function ROI_modulationAlgorithm_graph(roiHandle, selected_modality, selected_channel)

fitCell = roiHandle.UserData.fitRGB;
xCell = roiHandle.UserData.xCell;
yCell = roiHandle.UserData.yCell;
% legendCell = roiHandle.UserData.legendCell;

% Plot 1: Model, peak finding algorithm, and raw spectra

fitresult = fitCell{selected_modality}{selected_channel}; %modalities
xData = xCell{selected_modality}{selected_channel};
yData = yCell{selected_modality}{selected_channel};
%     legend_str = legendCell{k};

[s] = sinDescriptor(fitresult);
sinDescriptorPlot(fitresult, s, xData, yData)

% hFig.WindowState = 'maximized'; 

end