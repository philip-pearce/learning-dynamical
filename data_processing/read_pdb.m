villin = pdbread('villin_native.pdb')
xyz(:,1) = [villin.Model.Atom.X]';
xyz(:,2) = [villin.Model.Atom.Y]';
xyz(:,3) = [villin.Model.Atom.Z]';

atom = {villin.Model.Atom.AtomName}';
carbon = cellfun(@(x) strcmp(x,'CA'),atom);
xyzNative = xyz(carbon==1,:);
