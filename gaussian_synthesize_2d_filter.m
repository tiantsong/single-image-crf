function kernel_cell = gaussian_synthesize_2d_filter(L,half_window_n,scale)

% L = fitting polynomial order
% half_window_n = half of the window (excluding the middle sample)

N = 2*half_window_n+1;
count = 1;

% note : cell{x_order+1,y_order+1}

kernel_cell = cell(L+1,L+1);

for order = 0:L
    for n = 0:order

        xn = order-n;
        yn = n;
              
        hx = gaussian_synthesize_1d_filter(xn,scale,half_window_n);
        hy = gaussian_synthesize_1d_filter(yn,scale,half_window_n);
        
        ker = hy(:)*(hx(:)');
        
        % adapt for x-y axis
        kernel_cell{xn+1,yn+1} = flipud(ker);
        
        count = count + 1;        
    end
end

