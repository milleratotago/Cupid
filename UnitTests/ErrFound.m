function tf = ErrFound(sFName)

fid = fopen(sFName, 'r');
c1 = fscanf(fid,'%c');
fclose(fid);

% contains(str,pattern,'IgnoreCase',true)

% Error messages always to be ignored:
c1 = strrep(c1,'Warning: MATLAB has disabled some advanced graphics rendering features','');
c1 = strrep(c1,'Warning: Not enough discrete values to do probit estimation','');

% Error messages sometimes arising in computation of Fisher information:
c1 = strrep(c1,'Warning: Nearly singular information matrix','');
c1 = strrep(c1,'Warning: Matrix is singular to working precision','');

% Numerical errors in integrations:
c1 = strrep(c1,'Exiting: Maximum number of function evaluations has been exceeded','');
c1 = strrep(c1,'Warning: Reached the limit on the maximum number of intervals in use','');
c1 = strrep(c1,'Approximate bound on error is','');
c1 = strrep(c1,'Warning: Minimum step size reached near','');

% The following message arises with discrete distributions having low variance
% because the same value is at both percentiles defining the boundaries.
c1 = strrep(c1,'Warning: Cannot perform PctBoundsEstTest for this distribution because LowerX >= UpperX','');

% Often arises with Rician in dMATLABc
c1 = strrep(c1,'[Warning: NCX2INV did not converge.]','');

tf = contains(c1,'warning','IgnoreCase',true) || contains(c1,'error','IgnoreCase',true) || contains(c1,'Failed by','IgnoreCase',true) || contains(c1,'exiting','IgnoreCase',true);

end
