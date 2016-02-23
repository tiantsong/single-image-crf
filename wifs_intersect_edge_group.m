function wifs_intersect_edge_group(fileinfo)

disp('wifs_intersect_edge_group');

% Gamma image with bilateral vieville

ylimit = [0,1.6,40];
ylimit_rq = [0,1.8,40];
xlimit = [0,1,80];



wdataname = fileinfo.dataname;

% ----------

[point_hist,x_mid,y_mid,RQcell,RQweight,hist2d_1] = wifs_get_intersect_histogram_12(wdataname,fileinfo.linegroupname);

save(fileinfo.lineintersectname,'x_mid','y_mid','point_hist','hist2d_1');
save(fileinfo.signaturedata,'RQcell','RQweight');
