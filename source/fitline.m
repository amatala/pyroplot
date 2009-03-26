function [estimates, model] = fitline(xdata, ydata)
start_point = rand(1,2);
f = figure;
plot(1./xdata, ydata);
hold on;
model = @lin;
estimates = fminsearch(model, start_point);
[sse, FittedCurve] = model(estimates);
plot(1./xdata, FittedCurve, 'g');
function [sse, FittedCurve] = lin(params)
        A = params(1);
        E = params(2);
        FittedCurve = -log(A)+E./(8.3145).*xdata;
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end