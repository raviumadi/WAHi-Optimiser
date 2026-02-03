function plotCSV3D(csvFile)
% plotCSV3D  Plot 3D points from a CSV file with unique colours
%
% CSV format:
%   x,y,z
%   ...
%
% Example:
%   plotCSV3D('8mics_double_tetrahedron.csv')

    % --- load CSV ---
    T = readtable(csvFile);

    % % Basic checks
    % if ~all(ismember({'x','y','z'}, lower(T.Properties.VariableNames)))
    %     error('CSV must contain columns named x, y, z');
    % end

    % Extract coordinates
    X = T{:,1};
    Y = T{:,2};
    Z = T{:,3};

    N = numel(X);

    % --- plot ---
    figure('Color','w');
    ax = axes; %#ok<LAXES>
    hold(ax,'on');
    grid(ax,'on');
    axis(ax,'equal');
    view(ax,3);

    cmap = lines(N);   % unique colour per point

    for i = 1:N
        scatter3(ax, X(i), Y(i), Z(i), 120, cmap(i,:), 'filled', ...
            'MarkerEdgeColor','k','LineWidth',1.2);
        text(X(i), Y(i), Z(i), sprintf('  %d',i), ...
            'FontWeight','bold','FontSize',10,'Color',cmap(i,:));
    end

    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title(sprintf('3D points from %s', csvFile), 'Interpreter','none');

    hold(ax,'off');
end