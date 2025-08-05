function naxis(ax)

    if ~exist('ax','var')
        ax = gca;
    end
    set(ax,'ticklen',[0.04, 0.07],'FontName','Arial','FontSize',12,'TickDir','out','LineWidth',1.0,'Box','off');
    set(ax,'Xcolor','k','YColor','k','ZColor','k');
end