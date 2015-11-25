function result = writeMatrix(filename, A, seperator)
    result = false;
    [M, N] = size(A);
    fileID = fopen(filename,'W');
    for i=1:M
        for j=1:N
            if (j < N)
                fprintf(fileID,'%.20f%s', A(i,j), seperator);
            else
                fprintf(fileID,'%.20f\n', A(i,j));
            end
        end
    end
    result = true;
    fclose(fileID);
end