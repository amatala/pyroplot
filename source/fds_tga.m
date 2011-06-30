function [Mass, Temp] = fds_tga(template_file, fdsexe, CHID, Chrom, LogScaling,Par,Indx,dT)
% FDS_TGA   Calls FDS for TGA application using template
%
% [Mass, Temp] = fds_tga(template_file, fdsexe, CHID, Chrom, LogScaling[, Param,Indx])
%   template_file   name of the fds input template (string)
%   fdsexe          name of the fds executable (string)
%   CHID            fds job name (string)
%   Chrom           values of the variables (vector)
%   LogScaling      toggle the logarithmic scaling (vector of 0's and 1's)
%   Param           values of the parameters (vector, optional)
%   Indx            integer index for the file naming
%   dT              data type: 0 for TGA, 1 for HRR and 2 for MLR
% 
% FDS_TGA reads the contents of the template_file and replaces every string
% VAR#., where # is an integer, with a value of the corresponding item of
% the input vector Chrom. Additionally, all the strings PAR# are replaced
% by the values of the Par vector. The output is written to a temporary
% file and used to run Fire Dynamics Simulator (FDS).
%
% The ouput Mass is the time and mass of the sample material (profile file output) 
% and Temp contains the time and temperatures of the reference surface and
% sample surface. 
% 
% Example:
% [M,T]=fds_tga('Test1\tga.fds','c:\fds5\fds5.exe','tga',[],0)
% This command runs FDS calling tga.fds as it is.
%
if (nargin<7), Indx = 0; end
if (nargin<6), Par = []; end

%d = 2; %1 for burning rate (6. column), 2 for surface density (5. column)

%
% Create input file
%
namestr = [num2str(round(1000*cputime)) '_' num2str(Indx)];
data_file = ['fdstga' namestr '.fds'];
err_file = ['fdstga' namestr '.err'];
fid_in = fopen(template_file,'rt');
fid_out = fopen(data_file,'wt');
% print out the variables and parameters to data file
for i = 1:length(Chrom)
   s = sprintf('%0.8g',Chrom(i));
   fprintf(fid_out,'%s',[s ' ']);
end
fprintf(fid_out,'\n');
for i = 1:length(Par)
   s = sprintf('%0.8g',Par(i));
   fprintf(fid_out,'%s',[s ' ']);
end
fprintf(fid_out,'\n');
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
%
% Run FDS
%
sysstr = [fdsexe ' ' data_file ' 2> ' err_file];
system(sysstr);
%
% Read results
%
if strcmp(CHID,'tga') 
% devc: (1) time (2) front temperature (3) back temperature (4) surface density 
massfile = [CHID '_devc.csv'];
M = readdata(massfile,2);
Mass = M(:, [1, 4]); % time, surface density
Temp = M(:,[1, 2, 3]);
%M = M(:,[1 4 5]);
%Mass = [M(:,1),M(:,2).*M(:,3)];
%tempfile = [CHID '_devc.csv'];
%Temp = readdata(tempfile, 3);
% N = length(Mass);
% 
% l=1;
%      if min(Temp1(:,1))>min(Mass(:,1))  %if first value smaller        
%             while min(Temp1(:,1))>Mass(l,1)
%                 l=l+1;
%             end
%      end
%      
%     Temp(:,1) = Mass(l:length(Mass(:,1)),1);
%     Temp(:,2) = interp1(Temp1(:,1), Temp1(:,2), Mass(l:length(Mass(:,1)),1));
%     Temp(:,3) = interp1(Temp1(:,1), Temp1(:,3), Mass(l:length(Mass(:,1)),1));
%     Mass = Mass(l:length(Mass(:,1)),:);
%     
%     Temp(:,1) = removeNaNM(Temp(:,1),2);
%     Temp(:,2) = removeNaNM(Temp(:,2),3);
%     Temp(:,3) = removeNaNM(Temp(:,3),3);

elseif strcmp(CHID,'cone')
    file = [CHID '_hrr.csv'];
    file_devc = [CHID '_devc.csv'];
    data = readdata(file,2);
    data_devc = readdata(file_devc,2);
    HRR = data(:,2);
    if dT == 1
    MLR = interp1(data_devc(:,1),data_devc(:,6),data(:,1));
    Mass = [HRR, MLR]; %time, MLR
    elseif dT == 2  % if surface density!!
    g = -gradient(data_devc(:,5),data_devc(:,1));
    MLR = interp1(data_devc(:,1),g,data(:,1));
    Mass = [HRR, MLR]; %HRR, MLR
    end
    Temp = data(:,[2 1 3]); %HRR, time, rad_loss
    %end
end

if (1)
eval(['delete ' data_file])
eval(['delete ' err_file])
eval(['delete ' CHID '*.csv'])
eval(['delete ' CHID '*.out'])
eval(['delete ' CHID '*.end'])
eval(['delete ' CHID '*.smv'])
eval(['delete ' CHID '*.bf'])
eval(['delete ' CHID '*.sf'])
eval(['delete ' CHID '*.sz'])
eval(['delete ' CHID '*.s3d'])
end