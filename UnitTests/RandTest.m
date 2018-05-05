function thisval=RandTest(varargin)
    obj.N = 15;
    obj.P = .3;
    if numel(varargin)==0
        varargin{1} = 1;
    end
    unirands = rand(varargin{:},obj.N);
    thisvar01 = unirands>obj.P;
    thissize = size(thisvar01);
    thisval = squeeze(sum(thisvar01,numel(thissize)));
end



