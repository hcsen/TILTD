function plotperturbed(best, perturbed, lowerBound, upperBound, doHop)
    % Aux function to help with determining a sensible sigma value.
    % Replicates peturbation logic, if changed in main script must be
    % changed here also.
    [steps, length] = size(best);

    % First value is magnatudes different, is this right?
    x = 2:length;

    % If bounds given.
    if ~isempty(lowerBound) && ~isempty(upperBound)
        fill([x, fliplr(x)], [upperBound(2:end), fliplr(lowerBound(2:end))], [0.9 0.9 1], 'EdgeColor', 'none');
    end
            
    hold on;

    cmap = colormap(sky(steps));
    for i = 1:steps
        plot(x, best(i, 2:end), '-', 'Color', cmap(i,:));
        plot(x, perturbed(i, 2:end), '.', 'Color', cmap(i,:));   
    end
    axis tight;
   

    legend('Range', 'Best', 'Perturbed');
    grid on;

end