function ROI_modulation_graph(roiHandle, sel_modes, sel_modality)
%Input:
%sel_modes: range []
%sel_modality: 'RL PPL', 'RL XPL', 'TL PPL', 'TL XPL', 'TL XPL-lambda'

n_dataPoints = roiHandle.UserData.n_dataPoints;
n_pixels = roiHandle.UserData.n_pixels;
fitCell = roiHandle.UserData.fitRGB;
xCell = roiHandle.UserData.xCell;
yCell = roiHandle.UserData.yCell;
% legendCell = roiHandle.UserData.legendCell;

%% Plot 1: Raw spectra

lineWidth = 2;
fontSize = 13;
pointFactor = 7;

hFig = figure(9);
hFig.Position = [100, 100, 800, 600];
% hFig.Position = [100, 100, 1200, 800];

clf('reset') %clearing the automatic plot

%Plot grid setup

% Option 1: automatic (col>rows)

% nK = length(sel_modes);
% nf = ceil(nK^0.5); 
% if nf^2 - nK >= nf
%     nrows = nf-1;
%     ncolumns = nf;
% else 
%     nrows = nf;
%     ncolumns = nf;
% end
% clear nf

%Option 2: manual
nrows = 2;
ncolumns = 1;

t = tiledlayout(nrows, ncolumns, TileSpacing= "tight");
% %Nh, Nw, gap, marg_h, marg_w

k = 0;

for sel = sel_modes %sel_modes; 
    k = k + 1;
    m = k;

     % %Optional: skip grid tile    
    % if k >= 3
    %     m = m + 1;        
    % end

    fitresult = fitCell{sel}; %modalities
    xData = xCell{sel};
    yData = yCell{sel};   
    % legend_str = legendCell{k};    

    nexttile(m)
    
    %Channels
    h1 = plot(fitresult{1}, xData{1}, yData{1}, '.red');
    hold on
    h2 = plot(fitresult{2}, xData{2}, yData{2}, '.green');
    h3 = plot(fitresult{3}, xData{3}, yData{3}, '.blue');
    h4 = plot(fitresult{4}, xData{4}, yData{4}, '.black');
  
    % h1.MarkerSize = 10;
    %periods
    xline([0, 90, 180, 270, 360], ...
        'LineWidth', 0.8*lineWidth, 'Color', 'black', 'LineStyle','-') %period    
    hold off  

    %Axes
    y_unfolded = [yData{1}; yData{2}; yData{3}; yData{4}];    
    extraMargin = 0.05*max(y_unfolded, [], 'all'); %10    
    %xtop = max(xData{1})+ 1;
    xtop = 360;
    
    ymin = min(y_unfolded, [], 'all') - extraMargin;
    ymax = max(y_unfolded, [], 'all') + extraMargin;
    xlim([0, xtop]) %assuming equal x-range
    ylim([ymin, ymax])
%     ylim([0, 255])

    grid on
    ax = gca;
    grid(ax,'minor')
    legend(ax, 'off')
    % legend(legend_str, 'Location', 'EastOutside', 'Interpreter', 'none');

    %Tunning appeareance
    set(h1(1), 'MarkerSize', pointFactor*lineWidth, 'DisplayName', 'Red')
    set(h2(1), 'MarkerSize', pointFactor*lineWidth, 'DisplayName', 'Green')
    set(h3(1), 'MarkerSize', pointFactor*lineWidth, 'DisplayName', 'Blue')
    set(h4(1), 'MarkerSize', pointFactor*lineWidth, 'DisplayName', 'Greyscale')
    set(h1(2), 'color', [1, 0, 0, 0.5], 'LineWidth', lineWidth, 'DisplayName', 'R fit')
    set(h2(2), 'color', [0, 1, 0, 0.5], 'LineWidth', lineWidth, 'DisplayName', 'G fit')
    set(h3(2), 'color', [0, 0, 1, 0.5], 'LineWidth', lineWidth, 'DisplayName', 'B fit')
    set(h4(2), 'color', [0, 0, 0, 0.5], 'LineWidth', 1.5*lineWidth, 'DisplayName', 'Grey fit')    

    xlabel('Angle°', 'Interpreter', 'none', 'FontWeight','bold');
    ylabel('Channel 8-bit', 'Interpreter', 'none', 'FontWeight','bold');  
    title(sel_modality{sel}, 'FontSize', fontSize) %optional

    ax.XAxis.FontSize = 1*fontSize;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontSize = 1*fontSize;    
    ax.YAxis.FontWeight = 'bold';
end

%Optional
% leg = legend([h1(1), h2(1), h3(1), h4(1), h1(2), h2(2), h3(2), h4(2)]);
% leg.Layout.Tile = 'east';
% leg.Title.String = {sprintf('Rotation points = %d', n_dataPoints), sprintf('ROI pixels = %d', n_pixels)};

%Hiding axes labels
[row, col] = tilerowcol(t.Children);
if (max(row) > 1)
    xticklabels(t.Children(row < t.GridSize(1)),"") 
    xlabel(t.Children(row < t.GridSize(1)),"")        
end
if (max(col) > 1)
    % yticklabels(t.Children(col > 1), "") %do not activate
    ylabel(t.Children(col > 1), "")
end

% hFig.WindowState = 'maximized'; 

end