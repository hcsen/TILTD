function plotperturbed(lowerBound, upperBound, best, sigmas, MBH_tail, MBH_theta, MBH_scale)
    % Aux function to help with determining a sensible sigma value.
    % Replicates peturbation logic, if changed in main script must be
    % changed here also.
    size = length(best);

    pm = fix(rand(1, size) + 0.5);
    pm(~pm) = -1;
    perturbed = best + gprnd(MBH_tail, MBH_scale*sigmas, MBH_theta) .* pm;

    % First value is magnatudes different, is this right?
    x = 2:size;
    fill([x, fliplr(x)], [upperBound(2:end), fliplr(lowerBound(2:end))], [0.9 0.9 1], 'EdgeColor', 'none');
    hold on;
    plot(x, best(2:end), '-');
    plot(x, sigmas(2:end), '.');
    axis tight;
 
    plot(x, perturbed(2:end), 'x');      

    legend('Range', 'Best', 'Sigmas', 'Perturbed');
    grid on;

end