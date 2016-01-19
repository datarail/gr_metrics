function plot_GR_results(t_GRvalues, t_GRmetrics)

cell_lines = unique(t_GRvalues.cell_line);
agents = unique(t_GRvalues.agent);
xc = 10.^(-3.5:.01:1.5);
colors = [.8 .1 0; .7 1 .1; .2 0 .8];
for ia = 1:length(agents)
    % one figure window per drug
    figure(ia);clf;
    set(gcf,'position', [ia*30 ia*30 800 700])
    
    agent = agents{ia};
    for ip = [0 1]
        for it = [48 72]
            for ir = 1:3
                subplot(3,4, (ir-1)*4+ip*2+(it-24)/24)
                title(sprintf('t=%ih, perturb=%i', it, ip))
                hold on
                h = [];
                for iC = 1:length(cell_lines)
                    
                    t_ = t_GRvalues(t_GRvalues.perturbation==ip & ...
                        t_GRvalues.time==it & t_GRvalues.replicate==ir & ...
                        strcmp(t_GRvalues.cell_line, cell_lines{iC}) & ...
                        strcmp(t_GRvalues.agent,agent),:);
                    if isempty(t_), continue, end
                    h(iC) = plot(log10(t_.concentration), t_.GRvalue, 'o', 'color', ...
                        colors(iC,:));
                    
                    
                    t_ = t_GRmetrics(t_GRmetrics.perturbation==ip & ...
                        t_GRmetrics.time==it & t_GRmetrics.replicate==ir & ...
                        t_GRmetrics.cell_line==cell_lines{iC} & ...
                        t_GRmetrics.agent==agent,:);
                    
                    plot(log10(xc), t_.GRinf + (1-t_.GRinf)./ ...
                        ( 1 + (xc/t_.EC50).^t_.Hill), '-', 'color', ...
                        colors(iC,:))
                end
                xlim([-3.5 1.5])
                ylim([-1 1.3])
                set(gca,'xtick', -3:1, 'xticklabel', 10.^(-3:1), 'fontsize', 7)
                if ip==0 && it==48 
                    ylabel('GRvalue')
                end
                if ir==3
                    xlabel('concentration')
                end
                
                if ((ir-1)*4+ip*2+(it-24)/24)==1
                    hl = legend(h, cell_lines, 'orientation', 'horizontal');
                    set(hl,'position',[.1 .01 .8 .04])
                end
            end
        end
    end
end