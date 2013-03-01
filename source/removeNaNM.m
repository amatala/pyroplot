function [x1, y1, x2, y2, y3] = removeNaNM(x1, y1, x2, y2, y3)

%xy 1 - exp
%xy 2 - model

%remove NaN from x1
f = find(~isnan(x1));
x1 = x1(f);
y1 = y1(f);

%interpolate values of y2 to x1
y2 = interp1(x2, y2, x1);
y3 = interp1(x2, y3, x1);

%remove NaN from all
f = find(~isnan(y2));
x1 = x1(f);
x2 = x1;
y1 = y1(f);
y2 = y2(f);
y3 = y3(f);

%check
if any(isnan(x1)) || any(isnan(x2)) || any(isnan(y1)) || any(isnan(y2))
   msgbox('All NaN not removed!') 
end

end