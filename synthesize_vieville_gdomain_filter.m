function kernel_cell = synthesize_vieville_gdomain_filter(polyOrder,half_window_n)

% kernel_cell = synthesize_vieville_gdomain_filter(polyOrder,half_window_n)

% polyOrder = fitting polynomial order
% half_window_n = half of the window (excluding the middle sample)


% ----------- spatial weight ------------

std_domain = half_window_n;
[X,Y] = meshgrid(-half_window_n:half_window_n,half_window_n:-1:-half_window_n);
weight_domain = exp(-(X(:).^2 + Y(:).^2)/(2*std_domain^2));

% ----------------------------------------

N = 2*half_window_n+1;

P = [];

% To follow the scan order of matrix for (:)
% [1 3]
% [2 4]

for x = -half_window_n:half_window_n
    for y = half_window_n:-1:-half_window_n
        
        row = [];

        % 1 st order: x y
        % 2 nd order: xx xy yy

        for order = 0:polyOrder
            for n = 0:order

                xn = order-n;
                yn = n;

                row = [row, pixel_coeff(xn,x)*pixel_coeff(yn,y)];

            end
        end

        P = [P ; row];
    end
end

weight_domain = weight_domain/sum(sum(weight_domain));
weight_mat = diag(weight_domain);

A = inv(P'*weight_mat*P)*P'*weight_mat;


% ---------- reorder kernel in to Cell ---------

count = 1;

% cell{x_order+1,y_order+1}

kernel_cell = cell(polyOrder+1,polyOrder+1);

for order = 0:polyOrder
    for n = 0:order

        xn = order-n;
        yn = n;
              
        ker = A(count,:);
        ker = reshape(ker,N,N);        
        
        kernel_cell{xn+1,yn+1} = ker;
        
        count = count + 1;
    end
end

    function pc = pixel_coeff(fk,fi)
        
        pc = ((fi+0.5)^(fk+1) - (fi-0.5)^(fk+1))/factorial(fk+1);
        
    end
end