function ff = removeNan_2(F)

[r,s] = size(F);
ff = [];
if s == 2
    t = F(:,1);
    m = F(:,2);
l = length(t);
ll = length(m);
f = [];
k = 1;
if l ~= ll
   msgbox('t and m must be of the same length');
   return
else
    for i = 1:l
        if isnan(t(i))
        elseif isnan(m(i))
            f(k) = [t(i), m(i-1)];
        else
         f(k,:) = [t(i), m(i)];
         k = k+1;
        end
            
     end
        
end
ff = f(:,2);
elseif s == 1
    f = [];
    k = 1;
    for i = 1:r
        if ~isnan(F(i))      
           f(k) = F(i);
           k = k+1;
        end
    end
 ff = f';
end
    
end