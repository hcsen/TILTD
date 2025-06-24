function [goodRuns] = filterruns(optim, mass, violation)
    % Filters out implausable soultions via the following criteria.
    
    VIOLATION_THRESHHOLD = 0.1;
    count = struct( 'violation', 0 );
    [runs, length] = size(optim);
    goodRuns = zeros(0,length);
    for i = 1:runs()
        if violation(i) > VIOLATION_THRESHHOLD
            count.violation = count.violation +1;
            continue
        end
        goodRuns =  [ goodRuns; optim(i,:) ];
    end
    % Drop empty.
end

