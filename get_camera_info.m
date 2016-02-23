function [curvedata_nov08,curvedata_apr09,curvedata_sk] = get_camera_info(camera_id,color_id)

root = '/Users/tiantsong/Experiments/Camera-response-function/color-checker';

cam_C = {'canon_g3','canon_rebelxt','nikon_d70','kodak_dc290','sony_dcsv1','canon_450d','nikon_d90'};
cam_name = cam_C{camera_id};
useStandardRes = true;
Qorder = 'linear';
Qtype = 'pairwise';

% ---------- curve nov 08

curveDir = fullfile(root,'curves'); %'C:\projects\Camera-response-function\color-checker\curves';
fit_order = 1;
macbethFile = fullfile(curveDir,sprintf('%s_curves_%d.mat',cam_name,fit_order));

if exist(macbethFile,'file')
    load(macbethFile,'g_curve_k');

    k = g_curve_k{color_id};

    [Q,R,r] = eval_inv_polynomial_exponent_function(k,useStandardRes);

    curvedata_nov08.R = R;
    curvedata_nov08.r = r;
    curvedata_nov08.k = k;

    curvedata_nov08.Q = exp_5_calibrate_Q(Q,R,Qtype,Qorder);
else
    curvedata_nov08 = [];
end

% ------------- curve apr 09

curveDir = fullfile(root,'curves_apr09'); %'C:\projects\Camera-response-function\color-checker\curves_apr09';

macbethFile = fullfile(curveDir,sprintf('good_%s_data.mat',cam_name));

if exist(macbethFile,'file')

    load(macbethFile,'g_curve_k');

    k = g_curve_k{color_id};

    [Q,R,r] = eval_inv_polynomial_exponent_function(k,useStandardRes);

    curvedata_apr09.R = R;
    curvedata_apr09.r = r;
    curvedata_apr09.k = k;

    curvedata_apr09.Q = exp_5_calibrate_Q(Q,R,Qtype,Qorder);
else
    curvedata_apr09 = [];
end

% ----------- curve sk

curveDir = fullfile(root,'sk_curves'); %'C:\projects\Camera-response-function\color-checker\sk_curves';
macbethFile = fullfile(curveDir,sprintf('good_%s_data_sk.mat',cam_name));

if exist(macbethFile,'file')

    load(macbethFile,'g_curve_k');

    curvedata_sk = eval_sk_curve_function(g_curve_k(color_id),useStandardRes);
    curvedata_sk.Q = exp_5_calibrate_Q(curvedata_sk.Q,curvedata_sk.R,Qtype,Qorder);
else
    curvedata_sk = [];
end










