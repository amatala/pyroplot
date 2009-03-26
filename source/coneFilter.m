function yy = coneFilter(t,y)
% t time
% y hrr/mlr

%first interpolate to 1 s intervals, then set filter with 15 s moving averages

t_min=min(t);
t_max=max(t);

if round(t_min)<t_min
   t_min = round(t_min)+1;
else
    t_min = round(t_min);
end

if round(t_max)>t_max
   t_max = round(t_max)-1;
else
    t_max = round(t_max);
end

t_new = [t_min:t_max];

y_new = interp1(t,y,t_new);

y_filter = filtNs(y_new',15);

yy = interp1(t_new,y_filter',t);

if isnan(yy(1))
   yy(1)=yy(2); 
end

if isnan(yy(length(yy)))
   yy(length(yy)) = yy(length(yy)-1);
end

end