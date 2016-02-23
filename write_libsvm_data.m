function write_libsvm_data(filename,labelvector,datacell)

fid = fopen(filename,'a');

for i = 1:length(labelvector)
    label = labelvector(i);
    data = datacell{i};
    dataLength = size(data,1);
    
    for j = 1:size(data,2)
        
        fprintf(fid,'%d ',label);
        datavector = data(:,j);
        
        for k = 1:dataLength
            fprintf(fid,'%d:%.8f ',k,datavector(k));
        end
        fprintf(fid,'\n');
    end
end

disp('closing ...');
ST = fclose(fid);