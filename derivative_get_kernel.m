function kernel_cell = derivative_get_kernel(name,L,half_window_n)

if ~exist('L','var')
    L = 4;
end

switch name
    case 'vieville'
        if ~exist('half_window_n','var')
            half_window_n = 4;
        end
        
        ver_info = ver('matlab');
        
        if str2num(ver_info.Version) >= 7
            kernel_cell = vieville_synthesize_filter_ver2(L,half_window_n);
        else
            kernel_cell = vieville_synthesize_filter_matlab6(L,half_window_n);
        end
        
    case 'vieville_gdomain'
        
        if ~exist('half_window_n','var')
            half_window_n = 4;
        end


        kernel_cell = synthesize_vieville_gdomain_filter(L,half_window_n);


    case 'meer weighted'
        half_window_n = 5;
        kernel_cell = meer_synthesize_2d_weighted_filter(L,half_window_n);

    case 'meer unweighted'
        half_window_n = 4;
        kernel_cell = meer_synthesize_2d_unweighted_filter(L,half_window_n);

    case 'gaussian'
        scale = 1;
        half_window_n = 4;
        kernel_cell = gaussian_synthesize_2d_filter(L,half_window_n,scale);

    case 'farid'
        Lf = 3;
        half_window_n = 4;
        kernel_cell = farid_synthesize_filter(Lf,half_window_n);
end


