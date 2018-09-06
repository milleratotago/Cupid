function StopFailed(failedTests)
    % Rerun failed tests with StopOnFailures
    
    import matlab.unittest.TestSuite  % avoids having to specify this as a qualifier on its method names:
    runner = matlab.unittest.TestRunner.withTextOutput();
    import matlab.unittest.plugins.StopOnFailuresPlugin
    runner.addPlugin(StopOnFailuresPlugin);
    runner.run(failedTests);
    
end

