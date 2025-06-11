function plotperturbed(lowerBound, upperBound, best, MBH_sigma, doHop)
    % Aux function to help with determining a sensible sigma value.
    % Replicates peturbation logic, if changed in main script must be
    % changed here also.
    size = length(best);

    range = upperBound - lowerBound;
    perturbed = normrnd(best, MBH_sigma .* range); 
    if doHop
        for j = 1:size
            perturbed(dtIndex(j)) = perturbed(dtIndex(j)) + 2 * dtHop * rand - dtHop; % Hop each phase tof
        end
    end

    fprintf("( %u / %u ) values generated outside range.\n", sum(perturbed < lowerBound | perturbed > upperBound), size)

    % First value is magnatudes different, is this right?
    x = 2:size;
    fill([x, fliplr(x)], [upperBound(2:end), fliplr(lowerBound(2:end))], [0.9 0.9 1], 'EdgeColor', 'none');
    hold on;
    plot(x, best(2:end), '-');
    % plot(x, sigmas(2:end), '.');
    axis tight;
 
    plot(x, perturbed(2:end), 'x');      

    legend('Range', 'Best', 'Perturbed');
    grid on;

end