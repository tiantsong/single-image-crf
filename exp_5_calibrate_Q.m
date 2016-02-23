function [Qc,dQc] = exp_5_calibrate_Q(Q,R,type,order)

switch order
    
    case 'constant'
        switch type
            case 'gammacali'
                c0 = 2.5;
        end
        C = c0;
    
    case 'linear'
        switch type
            case 'graddiff'
                c0 = 0.2242;
                c1 = 26.5956;

            case 'rmsediff'
                c0 = 0.5775;
                c1 = 13.0045;

            case 'rmsediff-k2'
                c0 = 0.6818;
                c1 = 3.3218;
                
            case 'pairwise'
                c0 = 2.69978;
                c1 = 0.6616;
                   
            case 'pairwise_1'
                % from e1_compute_calibration_pairwise_1
                c0 = 2.5614;
                c1 = 0.6143;
                
             case 'pairwise_2'
                % from e1_compute_calibration_pairwise_2
                c0 = 2.5930;
                c1 = 0.6486;
                                                   
        end
        C = c0 + c1*R;



    case 'quadratic'
        switch type
            case 'graddiff'
                c0 = 0.8444;
                c1 = 0.0010;
                c2 = 37.9895;

            case 'rmsediff'
                c0 = 0.8467;
                c1 = 4.5471;
                c2 = 12.5478;

            case 'rmsediff-k2'
                c0 = 0.9737;    
                c1 = -0.1780;
                c2 = 4.3547;
                
        end
        C = c0 + c1*R + c2*R.^2;

end

Qc = C.*Q./(C.*Q - Q + 1);

if nargout > 1
    dQc = (C.^2.*Q - C.*Q  + 1)./(C.*Q - Q + 1).^2;
end


%
% Qc = sqrt(3)/(sqrt(3)-1)*(1-sqrt(1./(2*Q+1)));
%
% % Qc = sqrt(3)/(sqrt(3)-1)*(Q./(1.5*Q+1));