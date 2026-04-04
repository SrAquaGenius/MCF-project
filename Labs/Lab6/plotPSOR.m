function u = plotPSOR(g, x, u, xCrit, yCrit)
    % Plot de g(x)
    plot(x, g(x), 'ro')
    hold on
    
    % Plot da solução u
    plot(x, u, 'b.')
    
    % Letras gregas
    greek = {'\alpha','\beta','\gamma','\delta','\epsilon','\zeta','\eta','\theta',...
             '\iota','\kappa','\lambda','\mu','\nu','\xi','\pi','\rho'};
    
    % Plot das linhas verticais
    yl = ylim;
    yLabel = 0 - 0.01*(yl(2) - yl(1));
    
    for k = 1:length(xCrit)
        plot([xCrit(k), xCrit(k)], [0, yCrit(k)], '--k')
        text(xCrit(k), yLabel, greek{mod(k-1,length(greek))+1}, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'Rotation', 0)
    end
    
    % Eixos e legenda
    xline(0, 'k')
    yline(0, 'k')
    legend('g - obstacle', 'u - solution')
    
    hold off
end