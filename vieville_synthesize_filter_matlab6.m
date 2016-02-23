function kernel_cell = vieville_synthesize_filter_matlab6(L,half_window_n)

% L = fitting polynomial order
% half_window_n = half of the window (excluding the middle sample)

N = 2*half_window_n+1;

mid = N/2;

P = [];

for x = 0:(N-1)
    for y = 0:(N-1)

        row = [];

        for order = 0:L
            for n = 0:order

                xn = order-n;
                yn = n;

                fk = xn;
                fi = x;
                coef_x = ((fi+1-mid)^(fk+1) - (fi-mid)^(fk+1))/factorial(fk+1);
                                
                fk = yn;
                fi = y;
                coef_y = ((fi+1-mid)^(fk+1) - (fi-mid)^(fk+1))/factorial(fk+1);
                
                row = [row, coef_x*coef_y];

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

%     function pc = pixel_coeff(fk,fi)
% 
%         pc = ((fi+1-mid)^(fk+1) - (fi-mid)^(fk+1))/factorial(fk+1);
% 
%     end
end