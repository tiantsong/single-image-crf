function [point_hist,x_mid,y_mid,RQcell,RQweight,hist2d_1] = get_intersect_histogram_12(wdataname,line_resname)

load(wdataname,'R','Q');
load(line_resname,'vI_number','total_line','compSize','profile_type_info');

% --------------- Remove abnormal line ----------------------

ylimit = [0 2 200];
xlimit = [0 1 100];

x_edge = linspace(xlimit(1),xlimit(2),xlimit(3));
y_edge = linspace(ylimit(1),ylimit(2),ylimit(3));

x_mid = (x_edge(1:end-1) + x_edge(2:end))/2;
y_mid = (y_edge(1:end-1) + y_edge(2:end))/2;

x_len = length(x_edge)-1;
y_len = length(y_edge)-1;

disp_step = 1000;

I = vI_number>0.6;
vI_vec = vI_number(I);
R_vec = R(I);
Q_vec = Q(I);

flat_line_N = 0;
increase_seq_N = 0;
sharp_raise_N = 0;
large_gap_N = 0;

line_cell = cell(total_line,1);

fprintf('total_line = %d\n',total_line);

% -------- first pass begin ---------------

disp('first round');

remove_set = false(total_line,1);

for i = 1:total_line
        
    if i >= disp_step
        fprintf('%d ',i);
        disp_step = disp_step + 1000;
        if rem(disp_step,10000) == 0
            fprintf('\n');
        end
    end    

    I = find(vI_vec == i);
    
    Rgroup = R_vec(I);
    Qgroup = Q_vec(I);
    
    vI_vec(I) = [];
    R_vec(I) = [];
    Q_vec(I) = [];

    [Rgroup,sI] = sort(Rgroup);
    Qgroup = Qgroup(sI);

    % remove increasing sequence
    
    if any(diff(Qgroup) > 0)
        remove_set(i) = true;
        increase_seq_N = increase_seq_N + 1;
        continue;
    end
    
    % remove sharp raises        
    
%     if any(abs(Rgroup(end)-Rgroup(1))>0.1)
%         remove_set(i) = true;    
%         sharp_raise_N = sharp_raise_N + 1;
%         continue;
%     end
    
    % remove large gap
    
    if any(abs(Rgroup(2:end) - Rgroup(1:end-1))>0.12)
        remove_set(i) = true;
        large_gap_N = large_gap_N + 1;
        continue;        
    end
    
    % remove flat line
    
    if any(abs((Qgroup(end)-Qgroup(1))/(Rgroup(end)-Rgroup(1)))<2)
        remove_set(i) = true;        
        flat_line_N = flat_line_N + 1;
        continue;
    end
        
    line_cell{i} = [Rgroup(:), Qgroup(:)];
end

fprintf('\n increase_seq_N = %d, flat_line_N = %d, large_gap_N = %d\n',...
    increase_seq_N, flat_line_N,large_gap_N);

line_cell = line_cell(~remove_set);
clear remove_set

total_line = length(line_cell);

fprintf('total_line = %d\n',total_line);

% --------------------- remove duplicate ------------------

RQcell = cell(1,total_line);
RQcorner = cell(1,total_line);
RQweight = ones(1,total_line);
duplicate_set = false(total_line,1);

for i = 1:total_line

    if length(I) < 2
        fprintf('length(I) < 2\n');
    end

    Rgroup = line_cell{i}(:,1);
    Qgroup = line_cell{i}(:,2);

    Ri = Rgroup(1):0.001:Rgroup(end);
    Ri = sort([Ri(:) ; Rgroup(:)]);
    Rlen = length(Rgroup);
    if  Rlen < 3
        Qi = interp1q(Rgroup(:),Qgroup(:),Ri(:));
    else
        polydata = polyfit(Rgroup,Qgroup,2);
        Qi = polyval(polydata,Ri);
    end

    % -------------

    [hist2d,x_mid,y_mid] = tools_hist2d (Ri, Qi, x_edge, y_edge);
    RQcell{i} = find(hist2d > 0);
    RQcorner{i} = [min(Ri) max(Ri) min(Qi) max(Qi)];
end

disp('RQcell computation done');

disp_step = 1000;

for i = 1:total_line

    if rem(i,1000) == 0
        fprintf('scan = %d(dup = %d) ',i,sum(duplicate_set));
        if rem(i,5000) == 0
            fprintf('\n');
        end
    end

    if duplicate_set(i)
        continue;
    end

    A = RQcell{i};
    Al = length(A);
    Ac = RQcorner{i};

    Jindex = (i+1):total_line;
    JindexI = duplicate_set(Jindex);
    Jindex = Jindex(~JindexI);
    
    for j = Jindex

        Bc = RQcorner{j};
        
        Intersect_SQ = max((Ac(2)-Ac(1))+(Bc(2)-Bc(1))-(max(Ac(2),Bc(2))-min(Ac(1),Bc(1))),0)*max((Ac(4)-Ac(3))+(Bc(4)-Bc(3))-(max(Ac(4),Bc(4))-min(Ac(3),Bc(3))),0);
        Union_SQ = (max(Ac(2),Bc(2))-min(Ac(1),Bc(1)))*(max(Ac(4),Bc(4))-min(Ac(3),Bc(3)));
        
        if Intersect_SQ/Union_SQ < 0.111 %1/9
            continue;
        end
        
        B = RQcell{j};
        Bl = length(B);
                
        C = intersect(A,B);
        Cl = length(C);
        
        U = union(A,B);
        Ul = length(U);

        if Cl/Ul > 0.4 %0.333 %1/3 
            RQcell{j} = C;
            RQcorner{j} = [max(Ac(1),Bc(1)) min(Ac(2),Bc(2)) max(Ac(3),Bc(3)) min(Ac(4),Bc(4))];
            RQweight(j) = max(RQweight(i),RQweight(j)) + (1 - Cl/Ul);
            duplicate_set(i) = true;
            break;
        end
    end
end

RQcell = RQcell(~duplicate_set);
RQweight = RQweight(~duplicate_set);
fprintf('#duplicate_set = %d\n',sum(duplicate_set));
clear duplicate_set RQcorner

total_line = length(RQcell);
fprintf('total_line = %d\n',total_line);


% ------------------ remove long line -----------

x_mid = x_mid(:)';
Xmat = x_mid(ones(size(y_mid)),:);

y_mid = y_mid(:);
Ymat = y_mid(:,ones(size(x_mid)));

longline_set = false(total_line,1);
line_mid = zeros(total_line,1);

sharp_raise_N = 0;

for i = 1:total_line

    C = RQcell{i};
    Rgroup = Xmat(C);
    Qgroup = Ymat(C);
    
    maxR = max(Rgroup);
    minR = min(Rgroup);

    maxQ = max(Qgroup);
    minQ = min(Qgroup);
    
    temp = (maxR + minR)/2;
    [temp, minI] = min(abs(x_mid - temp));
    line_mid(i) = minI;    
    
    if abs(maxR-minR)>0.2 | abs(maxQ-minQ)>1
        longline_set(i) = true;    
        sharp_raise_N = sharp_raise_N + 1;
    end
    
end

fprintf('\n sharp_raise_N = %d\n',sharp_raise_N);

RQcell = RQcell(~longline_set);
line_mid = line_mid(~longline_set);
RQweight = RQweight(~longline_set);

clear longline_set RQcorner

total_line = length(RQcell);
fprintf('total_line = %d\n',total_line);

% -------- resampling ---------------

disp('resampling');

[histN,BIN] = histc_weight(x_mid(line_mid),RQweight,x_edge);
ave_N = round(mean(histN(histN > 0)));

sI = find(histN > ave_N);

for i = sI
    
    tempI = find(line_mid == i);
    RQweight(tempI) = RQweight(tempI)*ave_N/histN(i);
    
end

% sI = find(histN > ave_N);
% for i = sI
%     
%     tempI = find(line_mid == i);
%     tempN = length(tempI);
%     rI = randperm(tempN);
%     cutN = ceil((histN(i)-ave_N)/histN(i)*tempN);
%     
%     RQweight_this = RQweight(tempI);
%     RQweight_this = RQweight_this(:);
%     [RQweight_this,ssI] = sort(RQweight_this,1,'descend');
%     tempI = tempI(ssI);
%     
%     tempI = tempI(rI(1:cutN));
%     RQcell(tempI) = [];
%     RQweight(tempI) = [];
%     line_mid(tempI) = [];
% end

total_line = length(RQcell);
fprintf('total_line = %d\n',total_line);

% ------------------- intersect ------------------

point_hist = zeros(length(y_mid),length(x_mid));

disp('intersect');

for i = 1:total_line

    sI = RQcell{i};
    point_hist(sI) = point_hist(sI) + RQweight(i);
end

% ---------------------------
% remove non-interacting point

N_interact = 2;
I = point_hist < N_interact;
point_hist(I) = 0;

% ------- add single points -------

I = find(vI_number == 0.5);
fprintf('0.5 point = %d\n',length(I));
R_vec = R(I);
Q_vec = Q(I);

[hist2d_1,x_mid,y_mid] = tools_hist2d (R_vec, Q_vec, x_edge, y_edge);




