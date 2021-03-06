function kernel_cell = meer_synthesize_2d_unweighted_filter(L,half_window_n)

% L = fitting polynomial order
% half_window_n = half of the window (excluding the middle sample)

N = 2*half_window_n+1;

kernel_cell = cell(L+1,L+1);

count = 1;

for order = 0:L
    for n = 0:order

        xn = order-n;
        yn = n;

        hx = meer_synthesize_unweighted_filter(xn,L,half_window_n);
        hy = meer_synthesize_unweighted_filter(yn,L,half_window_n);
        
        ker = hy(:)*(hx(:)');
        
        name = 'I';
        for i = 1:xn
            name = [name,'x'];
        end

        for i = 1:yn
            name = [name,'y'];
        end

        % adapt for x-y axis
        kernel_cell{xn+1,yn+1} = flipud(ker);
               
        count = count + 1;        
    end
end


