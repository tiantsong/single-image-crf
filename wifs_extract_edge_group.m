function wifs_extract_edge_group(fileinfo)

disp('wifs_extract_edge_group');

% Gamma image with bilateral vieville

ylimit = [0,1.6,40];
ylimit_rq = [0,1.8,40];
xlimit = [0,1,80];

exdataname = fileinfo.exdataname;

wdataname = fileinfo.dataname;


load(wdataname,'I');


[vI_number,compSize,profile_type_info] = wifs_get_line_groups_consistent_fast_8dir(exdataname,I);

total_line = length(compSize);
if any(compSize<2)
    fprintf('some compSize<2\n');
end

save(fileinfo.linegroupname,'vI_number','total_line','compSize','profile_type_info');
