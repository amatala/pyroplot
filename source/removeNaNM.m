function f = removeNaNM(f,N)

% N = 1 -> Mass
% N = 2 -> Time
% N = 3 -> Temperature

if N == 1
    for i=1:length(f)
        if isnan(f(i))
            f(i) = f(i-1);
        end
    end
elseif N == 2
    dt = f(2)-f(1);
    for i=1:length(f)
        if isnan(f(i))
            f(i) = f(i-1)+dt;
        end
    end
elseif N == 3
    dT = f(2)-f(1);
    for i=1:length(f)
        if isnan(f(i))
            f(i) = f(i-1)+dT;
        end
    end
end
end