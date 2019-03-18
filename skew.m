function skewed = skew(mat)

% provides the skew symmetric form of a matrix
%
    skewed = [0 -mat(3) mat(2); mat(3) 0 -mat(1); -mat(2) mat(1) 0];


end