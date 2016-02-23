function [Io, Ix, Iy, Ixx, Iyy, Ixy, Ixxx, Iyyy, Ixxy, Ixyy] = derivatives_xy(R_image,name,outputOrder,polyOrder)

kernel_cell = derivative_get_kernel(name,polyOrder);

order_x = 0;
order_y = 0;
kernel_o = kernel_cell{1+order_x,1+order_y};

order_x = 1;
order_y = 0;
kernel_x = kernel_cell{1+order_x,1+order_y};

order_x = 0;
order_y = 1;
kernel_y = kernel_cell{1+order_x,1+order_y};

order_x = 2;
order_y = 0;
kernel_xx = kernel_cell{1+order_x,1+order_y};

order_x = 0;
order_y = 2;
kernel_yy = kernel_cell{1+order_x,1+order_y};

order_x = 1;
order_y = 1;
kernel_xy = kernel_cell{1+order_x,1+order_y};

order_x = 3;
order_y = 0;
kernel_xxx = kernel_cell{1+order_x,1+order_y};

order_x = 2;
order_y = 1;
kernel_xxy = kernel_cell{1+order_x,1+order_y};

order_x = 1;
order_y = 2;
kernel_xyy = kernel_cell{1+order_x,1+order_y};

order_x = 0;
order_y = 3;
kernel_yyy = kernel_cell{1+order_x,1+order_y};

Io = imfilter(R_image,kernel_o,'symmetric','same','corr');

if outputOrder >= 1
    Ix = imfilter(R_image,kernel_x,'symmetric','same','corr');
    Iy = imfilter(R_image,kernel_y,'symmetric','same','corr');
else
    Ix = 0;
    Iy = 0;
end

if outputOrder >= 2
    Ixx = imfilter(R_image,kernel_xx,'symmetric','same','corr');
    Iyy = imfilter(R_image,kernel_yy,'symmetric','same','corr');
    Ixy = imfilter(R_image,kernel_xy,'symmetric','same','corr');
else
    Ixx = 0;
    Iyy = 0;
    Ixy = 0;
end

if outputOrder >= 3
    Ixxx = imfilter(R_image,kernel_xxx,'symmetric','same','corr');
    Ixxy = imfilter(R_image,kernel_xxy,'symmetric','same','corr');
    Ixyy = imfilter(R_image,kernel_xyy,'symmetric','same','corr');
    Iyyy = imfilter(R_image,kernel_yyy,'symmetric','same','corr');
else
    Ixxx = 0;
    Ixxy = 0;
    Ixyy = 0;
    Iyyy = 0;
end
