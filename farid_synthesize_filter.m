function kernel_cell = farid_synthesize_filter(L,half_window_n)

% L = the highest order of the derivative kernel (max 2)
% half_window_n = half of the window (excluding the middle sample) (up to
% 4)

tap_number = 2*half_window_n+1;
count = 1;

% note : cell{x_order+1,y_order+1}

kernel_cell = cell(L+1,L+1);

for order = 0:L
    for n = 0:order

        xn = order-n;
        yn = n;
              
        hx = farid_synthesize_1d_filter(xn,half_window_n);
        hy = farid_synthesize_1d_filter(yn,half_window_n);
        
        ker = hy(:)*(hx(:)');        
        
        % adapt for x-y axis
        kernel_cell{xn+1,yn+1} = flipud(ker);
               
        count = count + 1;               
        
%         name = '';
% 
%         for i = 1:xn
%             name = [name,'x'];
%         end
% 
%         for i = 1:yn
%             name = [name,'y'];
%         end
%         
%         if xn+yn == 0
%             name = 'o';
%         end
%                
%         kernel_cell{xn+1,yn+1} = farid_generate_derivative_kernel(name,tap_number);
%         
%         count = count + 1;
    end
end

