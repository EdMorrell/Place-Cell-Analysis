function [] = RasterPlotter(spike_pos_iRun,timestamps,s_mean_frs,PC_Ver)

%   ---   RasterPlotter
%            - Function plots a raster for each cell (or PC)
%
% Inputs:
%        - spike_pos_iRun: Generated using get_SpikePos_IndRuns.m
%        - timestamps: Run timestamps in form [start end]
%        - s_mean_frs: Smoothed (can also be non-smoothed) mean firing
%                      rates for each cell
%        - PC_Ver: Array with 1 for every PC and 0 for non-PC's (default is
%                  plots all cells)

if nargin < 4
    PC_Ver = ones(size(spike_pos_iRun,1),1);
end

%%
for iCell = 1:size(spike_pos_iRun,1)
    if PC_Ver(iCell) == 1
        figure;
        subplot(2,1,1)
        hold on
        for iRun = 1:size(timestamps,1)
            if isempty(spike_pos_iRun{iCell,1}{iRun,1})
                continue
            else
                for iSpike = 1:size(spike_pos_iRun{iCell,1}{iRun,1},1)
                    plot([spike_pos_iRun{iCell,1}{iRun,1}(iSpike,1) ...
                    spike_pos_iRun{iCell,1}{iRun,1}(iSpike,1)],[(iRun-0.25)...
                    (iRun+0.25)],'c')
                    %Plot Properties
                    ax = plot_prop();
                end
            end
        end        
        ylim([0 (size(timestamps,1) + 1)])
        xlim([0 725])

        subplot(2,1,2)

        plot(s_mean_frs(iCell,:),'LineWidth',0.001)
        area(s_mean_frs(iCell,:),'FaceColor','c','EdgeAlpha',0)
        alpha(0.25)
        set(gca,'XTickLabel',[]);
        set(gca,'XTick',[])
        %Plot Properties
        ax = plot_prop();
        xlim([1 29])
        ylim([0 inf])

    end
end