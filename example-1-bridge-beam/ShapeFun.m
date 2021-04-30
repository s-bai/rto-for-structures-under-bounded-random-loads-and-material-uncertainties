function ShapeFunValue = ShapeFun(local_coordinate)
%% Shape function of 4-node bi-linear element
% Input:
%   local_coordinate: Local coordinates

% Output:
%   ShapeFunValue: Shape function values in row vector form

%   The shape functions are:
%   N1(xx, yy) = 1/4*(1 - xx)*(1 - yy);
%   N2(xx, yy) = 1/4*(1 + xx)*(1 - yy);
%   N3(xx, yy) = 1/4*(1 + xx)*(1 + yy);
%   N4(xx, yy) = 1/4*(1 - xx)*(1 + yy);

%%
xx = local_coordinate(1);
yy = local_coordinate(2);

% ShapeFunValue = [
%     1/4*(1 - local_coordinate(1))*(1 - local_coordinate(2)),...
%     1/4*(1 + local_coordinate(1))*(1 - local_coordinate(2)),...
%     1/4*(1 + local_coordinate(1))*(1 + local_coordinate(2)),...
%     1/4*(1 - local_coordinate(1))*(1 + local_coordinate(2))];

% ShapeFunValue_ori = [N1(xx, yy), N2(xx, yy), N3(xx, yy), N4(xx, yy)];


ShapeFunValue = 1/4*[(1 - xx)*(1 - yy), (1 + xx)*(1 - yy), (1 + xx)*(1 + yy), (1 - xx)*(1 + yy)];

end % End of ShapeFun()

%% Shape functions
% function N1_value = N1(xx, yy)
% %% N1
% N1_value = 1/4*(1 - xx)*(1 - yy);
% end
% 
% function N2_value = N2(xx, yy)
% %% N2
% N2_value = 1/4*(1 + xx)*(1 - yy);
% end
% 
% function N3_value = N3(xx, yy)
% %% N3
% N3_value = 1/4*(1 + xx)*(1 + yy);
% end
% 
% function N4_value = N4(xx, yy)
% %% N4
% N4_value = 1/4*(1 - xx)*(1 + yy);
% end

