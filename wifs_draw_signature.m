function wifs_draw_signature(fileinfo)

ylimit = [0,1.6,40];
ylimit_rq = [0,1.8,40];
xlimit = [0,1,80];

switch fileinfo.curvetype
    case 'camera'

        [curvedata_nov08,curvedata_apr09,curvedata_sk] = get_camera_info(fileinfo.camera_id,fileinfo.color_id);

        R = curvedata_sk.R;
        Q = curvedata_sk.Q;

    case 'gamma'

        load('B.mat','B');
        R = B;
        ga = exp_5_calibrate_Q(fileinfo.gamma_id,1,fileinfo.Qtype,fileinfo.Qorder);
        Q = ga*ones(size(R));

end


xEdge = linspace(xlimit(1),xlimit(2),xlimit(3));
yEdge = linspace(ylimit(1),ylimit(2),ylimit(3));

load(fileinfo.lineintersectname,'x_mid','y_mid','point_hist');

figure;
subplot(4,2,1); hold on;


imagesc(x_mid,y_mid,point_hist);
plot(R,Q,'-r','linewidth',2);
switch fileinfo.curvetype
    case 'camera'

        title(sprintf('R-Q: cam = %d',fileinfo.camera_id));
    case 'gamma'
        title(sprintf('R-Q: g = %.2f',fileinfo.gamma_id));
end

axis xy;
ylim([0 1.6]);
xlim([0 1]);
