function setup_time_selector(s, t, U)

idx_list = round(linspace(1, length(t), 5));

idx = idx_list(1);
hLine = plot3(s, t(idx)*ones(size(s)), U(idx,:), ...
              'k', 'LineWidth', 2);

for k = 1:length(idx_list)
    i = idx_list(k);
    
    uicontrol('Style','pushbutton',...
        'String',['t = ' sprintf('%.2f', t(i))],...
        'Position',[10 30*k 100 25],...
        'Callback', @(~,~) updateLine(i));
end

    function updateLine(i)
        set(hLine, ...
            'YData', t(i)*ones(size(s)), ...
            'ZData', U(i,:));
    end

end