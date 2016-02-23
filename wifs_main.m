function wifs_main

% Copyright: Tian-Tsong Ng, Jan 2009


cam_name_set = {'Canon G3','Canon Rebelxt','Nikon D70','camera4','Sony DCSv1'};

imdir = 'images';
datadir = 'data';

% set fileinfo.gamma_id and fileinfo.camera_id below according to image
% type

imname = 'gamma_1.0'; 
%imname = 'sony_dcsv1';

% image file path

fileinfo.imname = fullfile(imdir,[imname,'.tif']);

% R and Q data for selected points in an image

fileinfo.dataname = fullfile(datadir,[imname,'_RQdata.mat']);

% derivative data

fileinfo.exdataname = fullfile(datadir,[imname,'_diffdata.mat']);

% line group data

fileinfo.linegroupname = fullfile(datadir,[imname,'_linegroupname.mat']);

% intersection of line group

fileinfo.lineintersectname = fullfile(datadir,[imname,'_lineintersectname.mat']);

% signature data

fileinfo.signaturedata = fullfile(datadir,[imname,'_signaturedata.mat']);

% camera_id
% 1 = canon g3,
% 2 = canon rebelxt,
% 3 = nikon d70,
% 5 = sony dscv1,
% 0 = non-camera

fileinfo.camera_id = 0;

% gamma_id
% valid value range : 0 < gamma_id <= 1
% non-gamma image : 0

fileinfo.gamma_id = 1.0;

% color_id (data provided is for R channel)
% 1 = R
% 2 = G
% 3 = B

fileinfo.color_id = 1;

if fileinfo.camera_id > 0

    fileinfo.camname = cam_name_set{fileinfo.camera_id};
    fileinfo.curvetype = 'camera';
    fileinfo.Qorder = 'linear';
    fileinfo.Qtype = 'pairwise';

elseif fileinfo.gamma_id > 0

    fileinfo.curvetype = 'gamma';
    fileinfo.Qorder = 'constant';
    fileinfo.Qtype = 'gammacali';

end


% compute differential quantities

wifs_compute_differential_quant(fileinfo);

% extract edge groups

wifs_extract_edge_group(fileinfo);

% discard duplicated line groups and intersect the non-duplicated ones

wifs_intersect_edge_group(fileinfo);

% draw the histogram for the line group intersection points

wifs_draw_signature(fileinfo);