function filename = create_input(template_file, Chrom, LogScaling,Par)
%   template_file   name of the fds input template (string)
%   Chrom           values of the variables (vector)
%   LogScaling      toggle the logarithmic scaling (vector of 0's and 1's)
%   Param           values of the parameters (vector, optional)

filename = 'bestInput.fds';

fid_in = fopen(template_file,'rt');
fid_out = fopen(filename,'wt');

%
while (~feof(fid_in))
    s = fgets(fid_in);
    % insert Chrom variables
    for i = 1:length(Chrom)
        if (LogScaling(i))
            VarVal = 10^Chrom(i);
        else
            VarVal = Chrom(i);
        end
        VarStr = sprintf('%0.8g',VarVal);
        slen = length(s);
        key = ['VAR' num2str(i) '.'];
        keylen = length(key);
        keyindx = strfind(s,key);
        if (~isempty(keyindx))
            s = [s(1:keyindx-1) VarStr s(keyindx+keylen:slen)];
        end
    end
    % insert Par parameters
    for i = 1:length(Par)
        ParVal = Par(i);
        ParStr = sprintf('%0.8g',ParVal);
        slen = length(s);
        key = ['PAR' num2str(i)];
        keylen = length(key);
        keyindx = strfind(s,key);
        if (~isempty(keyindx))
            s = [s(1:keyindx-1) ParStr s(keyindx+keylen:slen)];
        end
    end
    fprintf(fid_out,'%s',s);
end
fclose(fid_in);
fclose(fid_out);

end