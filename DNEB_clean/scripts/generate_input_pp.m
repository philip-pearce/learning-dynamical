function generate_input_pp(gmm,varargin)
% This script takes the .mat file and makes it readable 
% for my fortran minimising code this should be ran in the DNEB_clean/data folder

%clear folder first
delete DNEB_clean/data/fortran/*.out
%also clear output folder
delete DNEB_clean/output/*.txt
obj = gmm;
if nargin == 2
objOriginal = varargin{1};
end

NComponent = size(obj.ComponentProportion,2);

weight = obj.ComponentProportion;
save DNEB_clean/data/fortran/weight.out weight -ascii -double;

meanvector = obj.mu;
save DNEB_clean/data/fortran/meanvector.out meanvector -ascii -double;

%out file for scaling symmetric Gaussians
    sig = [];
    for i = 1:NComponent
        A = obj.Sigma(:,:,i);
        sig = [sig, obj.Sigma(1,1,i)];
    
        dets(i) = det(2*pi*A);
        B = inv(A);
        e=['save DNEB_clean/data/fortran/matrix' num2str(i) '.out B -ascii -double'];
        eval(e);
    end
    
%out file for scaling non-spherical Gaussians
%if no original dim GMM supplied (i.e. no scaling necessary)
if nargin==1
        scaling = [];
    for i = 1:NComponent
        scaling = [scaling, 1];
    end
%if original dim GMM supplied (=> scaling necessary)
elseif nargin==2
        scaling = [];
    for i = 1:NComponent
        A = obj.Sigma(:,:,i);
        A2 = objOriginal.Sigma(:,:,i);
        scaling = [scaling, sqrt(det(2*pi*A))/sqrt(det(2*pi*A2))];
    end
end  
save DNEB_clean/data/fortran/sig.out sig -ascii -double;

save DNEB_clean/data/fortran/scalings.out scaling -ascii -double;

save DNEB_clean/data/fortran/determinants.out dets -ascii -double;
