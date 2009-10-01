function new_file = change_input_ambient(filename)

new_file = 'ox_limited_input.fds';

fid_in = fopen(filename,'rt');
fid_out = fopen(new_file,'wt');

misc_exists = false;

while (~feof(fid_in))
   s = fgets(fid_in);
   
   % check if MISC -line exists
   % if so, check if solid_phase_only = false and h_fixed = 0
   % if not, set them
   % assume that only one MISC-line exists
   
   key = '&MISC';
   keyindx = strfind(s,key);
   if (~isempty(keyindx))
        %s = [s(1:keyindx-1) VarStr s(keyindx+keylen:slen)];
        fprintf(fid_out,'%s','');
        %check if &MISC in many lines
        key_end = '/';
        end_indx = strfind(s,key_end);
        found = false;
        
        if (~isempty(end_indx))
            found = true;
        end
        while found == false
           s2 = fgets(fid_in); 
           fprintf(fid_out,'%s','');
           end_indx = strfind(s2,key_end);
           if (~isempty(end_indx))
               found = true;
           end
        end
        s = '&MISC SOLID_PHASE_ONLY = .TRUE., H_FIXED = 10 /';
        misc_exists = true;
   end
   
   if ~misc_exists %if not found yet
       key_matl = '&MATL'; %MISC should exist before Matl-lines
       matl_indx = strfind(s,key_matl);
       if (~isempty(matl_indx))
           s_misc = '&MISC SOLID_PHASE_ONLY = .TRUE., H_FIXED = 10 /';
           fprintf(fid_out,'%s\n',s_misc);
           fprintf(fid_out,'%s\n','');
           misc_exists = true;
       end
   end
   fprintf(fid_out,'%s',s);
end



end