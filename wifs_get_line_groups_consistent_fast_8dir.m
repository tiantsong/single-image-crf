function [vI_number,compSize,profile_type_info] = wifs_get_line_groups_consistent_fast_8dir(exdataname1,vI,Rx,Ry)

if nargin < 3
    load(exdataname1,'Rx','Ry');
end

line_count = 1;
vI_number = zeros(size(vI));

vI_loop = vI;

vI_temp = bwareaopen(vI_loop,2,8);
vI_number(vI_loop & ~vI_temp) = -1;
vI_loop = vI_temp;

direction_list = [180 135 90 45 0 -45 -90 -135 -180];

[I,J] = find(vI_loop);
direction_vec = 1000*ones(size(vI_loop));

for i = 1:length(I)
    
    Ii = I(i);
    Ji = J(i);        
    [dummy,minI] = min(abs(direction_list - atan2(Ry(Ii,Ji),Rx(Ii,Ji))*180/pi));
    direction_vec(Ii,Ji) = direction_list(minI);
end
    
clear Ry Rx

compN = 0;
compSize = [];

% ---------------

count = 1;
profile_type_info(count).direction = -45;
profile_type_info(count).index = [1 1];

[vI_number,compN,compSize] = process_m45(direction_vec == -45,vI_number,compN,compSize);

profile_type_info(count).index(2) = compN;
count = count + 1;

% ---------------

profile_type_info(count).direction = 135;
profile_type_info(count).index = [compN+1 compN+1];

[vI_number,compN,compSize] = process_m45(direction_vec == 135,vI_number,compN,compSize);

profile_type_info(count).index(2) = compN;
count = count + 1;

% ---------------


profile_type_info(count).direction = 45;
profile_type_info(count).index = [compN+1 compN+1];

[vI_number,compN,compSize] = process_45(direction_vec == 45,vI_number,compN,compSize);

profile_type_info(count).index(2) = compN;
count = count + 1;

% ---------------

profile_type_info(count).direction = -135;
profile_type_info(count).index = [compN+1 compN+1];

[vI_number,compN,compSize] = process_45(direction_vec == -135,vI_number,compN,compSize);

profile_type_info(count).index(2) = compN;
count = count + 1;

% ---------------

profile_type_info(count).direction = 0;
profile_type_info(count).index = [compN+1 compN+1];

[vI_number,compN,compSize] = process_0(direction_vec == 0,vI_number,compN,compSize);

profile_type_info(count).index(2) = compN;
count = count + 1;

% ---------------

profile_type_info(count).direction = 180;
profile_type_info(count).index = [compN+1 compN+1];

[vI_number,compN,compSize] = process_0(direction_vec == 180 | direction_vec == -180,vI_number,compN,compSize);

profile_type_info(count).index(2) = compN;
count = count + 1;

% ---------------

profile_type_info(count).direction = 90;
profile_type_info(count).index = [compN+1 compN+1];

[vI_number,compN,compSize] = process_90(direction_vec == 90,vI_number,compN,compSize);

profile_type_info(count).index(2) = compN;
count = count + 1;

% ---------------

profile_type_info(count).direction = -90;
profile_type_info(count).index = [compN+1 compN+1];

[vI_number,compN,compSize] = process_90(direction_vec == -90,vI_number,compN,compSize);

profile_type_info(count).index(2) = compN;
count = count + 1;

% ---------- scan diagonal -45 

function [vI_number,compN,compSize] = process_m45(bmap,vI_number,compN,compSize)

[rN,cN] = size(bmap);

vI_number_temp = zeros(size(vI_number));

for ci = 1:cN
    start_point = [1 ci];
    if cN-ci+1 <= rN
        end_point = [cN-ci+1 cN];
    else
        end_point = [rN ci+rN-1];
    end
    
    sI = sub2ind([rN,cN],start_point(1):end_point(1),start_point(2):end_point(2));
    scan = bmap(sI);
    scan = bwareaopen(scan,2,4);
    scanlabel = bwlabel(scan,4);    
    subcompN = max(scanlabel);

    K = scanlabel > 0;
    scanlabel = scanlabel(K);
    sI = sI(K);    
    scanlabel = scanlabel + compN;
        
    vI_number_temp(sI) = scanlabel;
    
    for i = compN+1:compN+subcompN
        compSize(i) = sum(scanlabel == i);
    end
        
    
    compN = compN + subcompN;
    clear scanlabel
end


for ri = 2:rN
    start_point = [ri 1];
    if rN-ri+1 <= cN
        end_point = [rN rN-ri+1];
    else
        end_point = [ri+cN-1 cN];
    end
    
    sI = sub2ind([rN,cN],start_point(1):end_point(1),start_point(2):end_point(2));
    scan = bmap(sI);
    scan = bwareaopen(scan,2,4);
    scanlabel = bwlabel(scan,4);    
    subcompN = max(scanlabel);

    K = scanlabel > 0;
    scanlabel = scanlabel(K);
    sI = sI(K);    
    scanlabel = scanlabel + compN;
        
    vI_number_temp(sI) = scanlabel;
    
    for i = compN+1:compN+subcompN
        compSize(i) = sum(scanlabel == i);
    end

    compN = compN + subcompN;
    clear scanlabel
end


K = vI_number_temp > 0;
vI_number(K) = vI_number_temp(K);

hm_pat = [0 1 0 ; 1 -1 1 ; 0 1 0];

K = bwhitmiss(K,hm_pat);
K = K & bmap;
vI_number(K) = 0.5;

clear K vI_number_temp


% ---------- scan diagonal 45 

function [vI_number,compN,compSize] = process_45(bmap,vI_number,compN,compSize)

[rN,cN] = size(bmap);

vI_number_temp = zeros(size(vI_number));


for ri = 1:rN
    start_point = [ri 1];
    if ri <= cN
        end_point = [1 ri];
    else
        end_point = [ri-cN+1 cN];
    end
    
    sI = sub2ind([rN,cN],start_point(1):-1:end_point(1),start_point(2):end_point(2));
    scan = bmap(sI);
    scan = bwareaopen(scan,2,4);
    scanlabel = bwlabel(scan,4);    
    subcompN = max(scanlabel);
    
    
    K = scanlabel > 0;
    scanlabel = scanlabel(K);
    sI = sI(K);    
    scanlabel = scanlabel + compN;
        
    vI_number_temp(sI) = scanlabel;
    
    for i = compN+1:compN+subcompN
        compSize(i) = sum(scanlabel == i);
    end
        
    
    compN = compN + subcompN;
    clear scanlabel
end


for ci = 2:cN
    start_point = [rN ci];
    if rN-cN+ci >= 1
        end_point = [rN-cN+ci cN];
    else
        end_point = [1 ci+rN-1];
    end
    
    sI = sub2ind([rN,cN],start_point(1):-1:end_point(1),start_point(2):end_point(2));
    scan = bmap(sI);
    scan = bwareaopen(scan,2,4);
    scanlabel = bwlabel(scan,4);    
    subcompN = max(scanlabel);

    K = scanlabel > 0;
    scanlabel = scanlabel(K);
    sI = sI(K);    
    scanlabel = scanlabel + compN;
        
    vI_number_temp(sI) = scanlabel;
    
    for i = compN+1:compN+subcompN
        compSize(i) = sum(scanlabel == i);
    end
    
    compN = compN + subcompN;
    clear scanlabel
end


K = vI_number_temp > 0;
vI_number(K) = vI_number_temp(K);

hm_pat = [0 1 0 ; 1 -1 1 ; 0 1 0];

K = bwhitmiss(K,hm_pat);
K = K & bmap;
vI_number(K) = 0.5;

clear K vI_number_temp


% ---------- scan diagonal 0 

function [vI_number,compN,compSize] = process_0(bmap,vI_number,compN,compSize)

[rN,cN] = size(bmap);

for ri = 1:rN
    start_point = [ri 1];
    end_point = [ri cN];
   
    sI = start_point(2):end_point(2);
    scan = bmap(start_point(1),sI);
    scan = bwareaopen(scan,2,4);
    scanlabel = bwlabel(scan,4);    
    subcompN = max(scanlabel);

    K = scanlabel > 0;
    scanlabel = scanlabel(K);
    sI = sI(K);    
    scanlabel = scanlabel + compN;
        
    vI_number(start_point(1),sI) = scanlabel;
    
    for i = compN+1:compN+subcompN
        compSize(i) = sum(scanlabel == i);
    end
        
    
    compN = compN + subcompN;
    clear scanlabel
end

% --------------------------------------------

function [vI_number,compN,compSize] = process_90(bmap,vI_number,compN,compSize)

[rN,cN] = size(bmap);

for ci = 1:cN
    start_point = [1 ci];
    end_point = [rN ci];

    sI = start_point(1):end_point(1);
    scan = bmap(sI,start_point(2));
    scan = bwareaopen(scan,2,4);
    scanlabel = bwlabel(scan,4);    
    subcompN = max(scanlabel);

    K = scanlabel > 0;
    scanlabel = scanlabel(K);
    sI = sI(K);    
    scanlabel = scanlabel + compN;
        
    vI_number(sI,start_point(2)) = scanlabel;
    
    for i = compN+1:compN+subcompN
        compSize(i) = sum(scanlabel == i);
    end
    
    compN = compN + subcompN;
    clear scanlabel
end

