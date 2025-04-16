
% Allowable time variance in percent.
RUNTIME_VARIANCE = 10;
OUTPUT_VARIANCE = 1;

% An array of cell arrays, where each cell array is a case to test.
% In the format {'name', 'config', 'violation', 'runtime (s)'}

cases = [{'basic',{'Inputs/verification_conf_base', 'hpc_conf'}, [ 0.21839974760232377937 ], 400},...
         {'mbh',{'Inputs/verification_conf_mbh', 'hpc_conf'}, [ 0.012120264285739645871, 0.063459270359921549076 ], 600}];


for i = 1:size(cases)
    tStart = tic;
    output = lofi_search(cases{2}{:});
    tEnd = toc(tStart);
    
    result = "FAIL";
    reason = 'Crashed';
    passOutput = false;
    referenceResults = cases{3};
    if ~isequal(size(referenceResults), size(output.violation_archive))
        reason = 'OUTPUT SHAPE MISMATCH';
    else
        variance = sum(referenceResults - output.violation_archive) / referenceResults; 
        if OUTPUT_VARIANCE < variance * 100
            reason = sprintf('Output results not within %%%u of reference case (%%%d).', OUTPUT_VARIANCE, variance*100)
        else
            passOutput = true;
        end
    end
    tReference = (cases{4} * ((RUNTIME_VARIANCE * 0.01)+1) )
    if tEnd > tReference
        reason = sprintf('Runtime exceeds %%%u of reference case (%us/%us).',  RUNTIME_VARIANCE, tEnd, tReference);
    elseif passOutput
        reason = "everything looks good";
        result = "PASS";
    end
    disp(sprintf('Case: %s - %s - %s\n', cases{1}, result, reason));
end
