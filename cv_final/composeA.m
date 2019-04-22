% function A = composeA(x1, x2)
% Compose matrix A, given matched points (X1,X2) from two images
% Input: 
%   -normalized points: X1 and X2 
% Output: 
%   -matrix A
function A = composeA(x1, x2)
    x1x = x1(1,:)';
    x1y = x1(2,:)';
    x2x = x2(1,:)';
    x2y = x2(2,:)';

    A = [x1x.*x2x x1x.*x2y x1x x1y.*x2x x1y.*x2y x1y x2x x2y ones(size(x1x))  ];


end
