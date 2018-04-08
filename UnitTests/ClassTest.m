function [tests, results] = ClassTest(utAny,Parallel,runner)

classtime = tic;

tests =  matlab.unittest.TestSuite.fromClass(utAny);

if Parallel
    results =  runInParallel(runner,tests);
else
    results =  run(tests);
end

classmins_elapsed = toc(classtime) / 60

end
