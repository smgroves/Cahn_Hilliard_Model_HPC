% function arrOut = laplacian(arrIn)
% % Calculates laplacian for a matrix using a 3x3 convolution with edge wrapping
% arrOut = -arrIn + ...
%     0.2(circshift(arrIn,1)+circshift(arrIn,-1)+circshift(arrIn,-1,2)+circshift(arrIn,1,2));
%     %  + ...
%     % 0.05*(circshift(arrIn,[1 1])+circshift(arrIn,[1 -1])+circshift(arrIn,[-1 1])+circshift(arrIn,[-1 -1]));
% end

function arrOut = laplacian(arrIn)
    [nx,ny] = size(arrIn);

    arrOut = zeros(nx,ny); %Initialize Laplacian
    for i = 1:nx
        for j = 1:ny
            if i > 1
                dadx_L = arrIn(i,j)-arrIn(i-1,j);
            else
                dadx_L = 0;
            end
            if i < nx
                dadx_R = arrIn(i+1,j)-arrIn(i,j);
            else
                dadx_R = 0;
            end
            if j > 1
                dady_B = arrIn(i,j)-arrIn(i,j-1);
            else
                dady_B = 0;
            end
            if j < ny
                dady_T = arrIn(i,j+1)-arrIn(i,j);
            else
                dady_T = 0;
            end
            arrOut(i,j) = (dadx_R-dadx_L + dady_T-dady_B);
        end
    end
    
    end