oneDist = Normal(0,1);
data = oneDist.Random(1000,1);
[mean(data), std(data,1)]

disp('Example of basic ML estimation')
EstFn = @oneDist.EstML;
fnParms = {data};
%                   mu   sigma
starting_points = [-0.1, 0.9; ...
                    0.1, 1.1];
[s,EndingVals,fval,exitflag,output,allstarts]= EstManyStarts(oneDist,EstFn,fnParms,starting_points);
s

disp('Example of MLcensored estimation')
Bounds = [-2, 2];   % define truncation bounds
truncData = data(data>=Bounds(1) & data<=Bounds(2));       % censor the data
nsTooExtreme = [sum(data<Bounds(1)), sum(data>Bounds(2))]  % count too-small and too-large observations
EstFn = @oneDist.EstMLcensored;
fnParms = [{truncData}, {Bounds}, {nsTooExtreme}];
%                   mu   sigma
starting_points = [-0.1, 0.9; ...
                    0.1, 1.1; ...
                    0.0, 1.0];
[s,EndingVals,fval,exitflag,output]= EstManyStarts(oneDist,EstFn,fnParms,starting_points);
s

disp('Example of MLcensored estimation with sigma constrained')
% For this example, suppose we want to consider only parameter combinations
% with sigma exactly 0.9, 1.0, or 1.1, and we want the best of those.
Bounds = [-2, 2];   % define truncation bounds
truncData = data(data>=Bounds(1) & data<=Bounds(2));       % censor the data
nsTooExtreme = [sum(data<Bounds(1)), sum(data>Bounds(2))]  % count too-small and too-large observations
EstFn = @oneDist.EstMLcensored;
parmCodes = 'rf';  % indicate that sigma is constrained
fnParms = [{truncData}, {Bounds}, {nsTooExtreme}, parmCodes];  % pass optional parmCodes parameter
%                   mu   sigma
starting_points = [ 0.0, 0.9; ...
                    0.0, 1.0; ...
                    0.0, 1.1];
[s,EndingVals,fval,exitflag,output]= EstManyStarts(oneDist,EstFn,fnParms,starting_points);
s


