
% Allowable time variance in percent.
RUNTIME_VARIANCE = 10;
OUTPUT_VARIANCE = 1;

% An array of cell arrays, where each cell array is a case to test.
% In the format {'name', 'config', 'violation', 'runtime (s)'}

cases = {{'basic',{'Inputs/verification_conf_base', 'hpc_conf'}, [ 0.21839974760232377937 ], 540},...
         {'mbh',{'Inputs/verification_conf_mbh', 'hpc_conf'}, [ 0.012120264285739645871; 0.063459270359921549076 ], 640}};

parfor i = 1:length(cases)
    [results(i), times(i)] = testcase(cases{i});
end


for i = 1:length(cases)
    result = "FAIL";
    reason = 'Crashed';
    passOutput = false;
    referenceResults = cases{i}{3};
    testResults = results(i).violation_archive;
    if ~isequal(size(referenceResults), size(testResults))
        reason = sprintf('OUTPUT SHAPE MISMATCH');
    else
        variance = sum(referenceResults - testResults) / referenceResults; 
        if OUTPUT_VARIANCE < variance * 100
            reason = sprintf('Output results not within %u%% of reference case (%u%%).', OUTPUT_VARIANCE, variance*100);
        else
            passOutput = true;
        end
    end
    tReference = (cases{i}{4} * ((RUNTIME_VARIANCE * 0.01)+1) );
    if times(i) > tReference
        reason = sprintf('Runtime exceeds %u%% of reference case (%us/%us).',  RUNTIME_VARIANCE, times(i), tReference);
    elseif passOutput
        reason = "Everything looks good";
        result = "PASS";
    end
    fprintf('%s - %s - %s\n', cases{i}{1}, result, reason);
    failedResults =  results(i);
    save('failed_output.mat', 'failedResults');
    fprintf('Wrote failed results to ''failed_output.mat''');
end 

function [result, time] = testcase(input)
    tStart = tic;
    result = lofi_search(input{2}{:});
    time = toc(tStart);
end
