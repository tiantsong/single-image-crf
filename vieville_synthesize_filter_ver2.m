function kernel_cell = vieville_synthesize_filter_ver2(L,half_window_n)

% ver2 : According to the pixel index (not the boundary of a square pixel),
% the central pixel is (0 0)

% ver2 should be equal to the ver1.

% L = fitting polynomial order
% half_window_n = half of the window (excluding the middle sample)

N = 2*half_window_n+1;

P = [];

for x = -half_window_n:half_window_n
    for y = -half_window_n:half_window_n

        row = [];

        for order = 0:L
            for n = 0:order

                xn = order-n;
                yn = n;

                row = [row, pixel_coeff(xn,x)*pixel_coeff(yn,y)];

            end
        end

        P = [P ; row];
    end
end

A = inv(P'*P)*P';

count = 1;

% cell{x_order+1,y_order+1}

kernel_cell = cell(L+1,L+1);

for order = 0:L
    for n = 0:order

        xn = order-n;
        yn = n;
              
        ker = A(count,:);
        ker = reshape(ker,N,N);        
        
        % adapt for x-y axis
        kernel_cell{xn+1,yn+1} = flipud(ker);
        
        count = count + 1;
    end
end

    function pc = pixel_coeff(fk,fi)
        
        pc = ((fi+0.5)^(fk+1) - (fi-0.5)^(fk+1))/factorial(fk+1);
        
    end
end