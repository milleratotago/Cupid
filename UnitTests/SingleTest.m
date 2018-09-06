function failedTests = SingleTest(sThisClass)
    
    global WantPlots
    global GlobalSkipAllEst
    
    %% **** START OF OPTIONS ****
    if numel(WantPlots) == 0
        WantPlots = true;
        % WantPlots = false;
    end
    
    if numel(GlobalSkipAllEst)==0
        GlobalSkipAllEst = false;
        % GlobalSkipAllEst = true;
    end
    
    Parallel = false;
    % Parallel = true;
    
    warning off backtrace;   % Display 1-line warning messages.
    % warning on backtrace;    % Display full warning messages.
    
    % **** END OF OPTIONS ****
    
    %% Nothing changes from here.
    
    % Some standard lines:
    import matlab.unittest.TestSuite  % avoids having to specify this as a qualifier on its method names:
    runner = matlab.unittest.TestRunner.withTextOutput();
    
    ThisClass = eval(['?ut' sThisClass]);       % Specify the class of unit tests to run.
    DiaryName = ['00Results/' sThisClass '.err'];  % Specify a name for the output diary file.
    if exist(DiaryName,'file')>0
        delete(DiaryName)  % Delete old version of err file
    end
    diary(DiaryName);
    
    time0 = tic;
    
    [tests, results] = ClassTest(ThisClass,Parallel,runner);
    
    % Note: When running parallel, MATLAB reports xxx seconds testing time, but it
    % appears to sum the times for different processors working in parallel.
    totalmins_elapsed = toc(time0) / 60
    
    s=evalc('disp(results)');
    wantloc=strfind(s,'Totals:');
    fprintf(['Summary of test result totals:\n' s(wantloc+8:end)]);
    
    diary('off');
    
    if ErrFound(DiaryName)
        CleanHTML(DiaryName);
        disp(['Errors in ' DiaryName ', so save it.']);
        
        % Rerun failed tests with StopOnFailures
        % Save & return these failed tests.
%         import matlab.unittest.plugins.StopOnFailuresPlugin
%         runner.addPlugin(StopOnFailuresPlugin);
        failedTests = tests([results.Failed]);
%         result2 = runner.run(failedTests);

    else
        disp(['No errors in ' DiaryName ', so remove it.']);
        delete(DiaryName);
        failedTests = [];
    end
    
end % SingleTest
