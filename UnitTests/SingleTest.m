%% Modify this script inside the **** OPTIONS **** section.

% Some standard lines:
import matlab.unittest.TestSuite  % avoids having to specify this as a qualifier on its method names:
runner = matlab.unittest.TestRunner.withTextOutput();
global WantPlots
global GlobalSkipAllEst


% **************** OPTIONS START HERE ****************
% Note: Further speed options are available in utGeneric.m
% at the points marked *** CONTROL SPEED HERE ***

% Choose whether you want some standard plots to be generated.
  WantPlots = true;
% WantPlots = false;

% Choose whether you want to skip the parameter estimation part of the unit testing
% (which may be rather slow, and which you may not care about if you are not going
% to do any parameter estimation).
  GlobalSkipAllEst = false;
% GlobalSkipAllEst = true;

% Choose whether you want unit tests to be run in parallel.
  Parallel = false;
% Parallel = true;

% Choose the verbosity of the warning messages when there are problems.
  warning off backtrace;   % Display 1-line warning messages.
% warning on backtrace;    % Display full warning messages.

sThisClass = 'StudRng';
ThisClass = eval(['?ut' sThisClass]);       % Specify the class of unit tests to run.
DiaryName = [sThisClass '.err'];  % Specify a name for the output diary file.
if exist(DiaryName,'file')>0
    delete(DiaryName)  % Delete old version of err file
end

% **************** OPTIONS END HERE ****************

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

