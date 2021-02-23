function M=readComplexMatrix(fileID)
    fid = fopen(fileID,'r');
    A = textscan(fid,'%s');
    fclose(fid);
    tempFile=['Complex2NormalForm_' fileID];
    fid = fopen(tempFile,'w');
    for i=1:501
        s=replace(A{1}{i},'+-','-');
        fprintf(fid,s);
        fprintf(fid,'\n');
    end
    fclose(fid);
    M=readmatrix(tempFile);
    delete Complex2NormalForm_**;
end