function [bestChrom,oV] = fdsga(estimates,data,template, FdsExe, weights, Aindex, limits, LogScaling, par, parInd, moisture)

% Parameters
%
% fdsgaA(data,template, FdsExe, weights, wA, var, limits, LogScaling, par)
% estimates = initial guess
% template = fds input
% FdsExe = fds exe file
% weights = weights of fitness function, default [985/1000, 15/1000]
% Aindex is vector with the indexes of variables A
% limits = minimum and maximum values for variables;
%     [min;max]    
%     in order A's, E's, residues, others
% LogScaling = 1 if variable is logarithm, 0 if not
% par = vector where 1. cell is 1 if TRamp is parameter 
%                    2. Temperature area (Final temperature - starting
%                    temperature)
%                    3. is 1 if TFin is parameter
%                    4. TFin
%                    5. -> index of residue variable 
% parInd = vector that tells the index of parameters


% TGA data format
%
% data(i).Type     'TGA'
% data(i).Rate     Heating rate [C/min]
% data(i).T        Temperature [C]
% data(i).Time     Time (s)
% data(i).M        Sample mass [kg]
% data(i).dMdt     Gradient of mass fraction [kg/s]

% cone data format
% data(1).Type     'Cone'
% data(1).T        Time (s)
% data(1).M        MLR [kg/m^2] or HRR [kW/m^2]
% data(1).A        Area of sample
% data(1).dataType 1 if HRR, 2 if MLR
display('GA is running');
msgbox('GA is running');

% if isempty(weights)
%     weights = [985/1000 15/1000];
% end

Ndata = length(data);

% check for forgotten stop file
if (exist('ga.stop','file'))
   msgbox('Warning: stop file exists');
end

%make sure TGA is mass fraction
if strcmp(data(1).Type,'TGA')
    for i=1:Ndata    
    data(i).M = data(i).M./max(data(i).M);
    end
end
%figures
bestIndividual = figure('Visible', 'off', ...
                    'Name', 'Best Individual', ...
                    'NumberTitle', 'off');
bestIndividual_2 = figure('Visible', 'off', ...
                    'Name', 'Best Individual', ...
                    'NumberTitle', 'off');
fitness = figure('Visible', 'off',...
                    'Name', 'OBJECTIVE VALUES',...
                    'NumberTitle', 'off');
                
warning('off','MATLAB:dispatcher:InexactMatch');
%create initial population
%Nind is the number of variables of each individual
%FieldDR is 2xNind vector of minimum and maximum values
%Parameters
NIND = 20; %!!!!!!!!!!!!!!later 20
FieldDR = limits;
iLogScaling = find(LogScaling);
FieldDR(:,iLogScaling) = log10(FieldDR(:,iLogScaling));
[mF, NVAR] = size(FieldDR); % size of FieldDR, NVAR is number of variables
% Define GA Parameters
GGAP = 0.8; % Generation gap 0,8
XOVR = 0.7; % Crossover rate 0,7
MUTR = 1/NVAR; % Mutation rate
MAXGEN = 1200; % Maximum no. of generations
INSR = 0.9; % Insertion rate
SUBPOP = 4; % No. of subpopulations
MIGR = 0.2; % Migration rate
MIGGEN = 20; % No. of gens / migration

%create initial population
%if estimates exist, add them to initial population
if isempty(estimates)
    Chrom = crtrp(SUBPOP*NIND, FieldDR);
else
Chrom = crtrp(SUBPOP*NIND-1, FieldDR);
Chrom(SUBPOP*NIND,:)=estimates;
end
objVal = objF(Ndata, data, Chrom, LogScaling, template, FdsExe, weights, Aindex, limits, par, parInd, moisture);
[bestEver,line] = min(objVal);
gen = 0;
set(fitness, 'Visible', 'on');
xlabel('Generation');
ylabel('Fitness value');
while gen < MAXGEN & ~exist('ga.stop','file') 
    %get Fitness Values
FitnV = ranking(objVal,[2,1],SUBPOP);
%get selection
SelCh = select('sus', Chrom, FitnV, GGAP, SUBPOP);
%then recombine
SelCh = recombin('recdis', SelCh, XOVR, SUBPOP); %xovsh - 0.7
%mutate offsprings
SelCh = mutate('mutbga',SelCh,FieldDR,MUTR, SUBPOP);
% Calculate objective function for offsprings
ObjVOff = objF(Ndata, data, SelCh, LogScaling, template, FdsExe, weights, Aindex, limits, par, parInd, moisture);
% Insert best offspring replacing worst parents
[Chrom, objVal] = reins(Chrom, SelCh, SUBPOP,[1, INSR], objVal, ObjVOff);

% Migrate individuals between subpopulations
if (rem(gen,MIGGEN) == 0)
[Chrom, objVal] = ...
migrate(Chrom, SUBPOP, [MIGR, 1, 1], objVal);
end

%plot best and average objective values
if gen == 0
[best(1), line] = min(objVal);
else
[best(gen), line] = min(objVal);
end

averagefitness = mean(objVal);
set(0, 'CurrentFigure', fitness);
if gen == 0
h=plot(gen, best(1), 'go');
else
h=plot(gen, best(gen), 'go');
end
hold on;
plot(gen, averagefitness, 'o');
hold on;

if ishandle(bestIndividual)
    delete(bestIndividual);
    bestIndividual = figure('Visible', 'off', ...
                    'Name', 'Best Individual', ...
                    'NumberTitle', 'off');
end

if ishandle(bestIndividual_2)
    delete(bestIndividual_2);
    bestIndividual_2 = figure('Visible', 'off', ...
                    'Name', 'Best Individual', ...
                    'NumberTitle', 'off');
end


for k=1:Ndata
    %if strcmp(data(1).Type,'TGA')
    if strcmp(data(1).Type,'TGA')
    parameters = zeros(1,length(parInd));
    index = 5;
    NR = par(5);
    L = length(par);
    for n=1:length(parInd)
    
     if parInd(n)==1
         %if tramp T parameter
            if par(1)==1
            parameters(n) = par(2)/data(k).Rate*60;
            %if residue parameter
            %if mass fraction parameter
            elseif parInd(n)==length(parInd)
             parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L))); 
             index = index+1;
             %if residue parameter
            else
            index = index+1;
            parameters(n)=1-Chrom(i,par(index));
            end
        elseif parInd(n)==2
            if par(1)==1 && par(3)==1
             parameters(n)= par(4)/data(k).Rate*60;
            elseif parInd(n)==length(parInd)
             parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L))); 
             index = index+1;
             %if residue parameter
            else
            index = index+1;
            parameters(n)=1-Chrom(i,par(index));
            end
        elseif parInd(n)==length(parInd) && parInd(n) > (2+par(5))
                parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L)));
                index = index +1;
        else
            index = index+1;
            parameters(n)=1-Chrom(i,par(index));
            %if mass fraction parameter
     
     end %if
    end%for n

    [M,T]=fds_tga(template,FdsExe,'tga',Chrom(line,:),LogScaling,parameters,[],0);
    M(:,2)=M(:,2)./max(M(:,2)); 
    dt = T(2,1)-T(1,1); %time step
    T0 = T(1,2);
        if max(T(:,2))<max(data(k).T)
            i=length(T(:,2));
            while T(i-1,2)<max(data(k).T)
                M(i,1)=M(i-1,1)+dt;
                M(i,2)=M(i-1,2);
                T(i,1)=T(i-1,1)+dt;
                T(i,2)=(T(i,1)-dt)*data(k).Rate/60+T0;
                i=i+1;
            end
            
        end
elseif strcmp(data(1).Type,'Cone')
    
    if strcmp(data(k).gas,'Air')
        temp = template;
    else %if in N2
        temp = 'ox_limited_input.fds';
    end
    parameters = [];
    [M,T]=fds_tga(temp,FdsExe,'cone',Chrom(line,:),LogScaling,parameters,[],data(k).dataType);
    M(:,1)=coneFilter(T(:,2),M(:,1));
    M(:,2)=coneFilter(T(:,2),M(:,2));
    M(:,1)=M(:,1)./data(1).A;
end

set(0,'CurrentFigure', bestIndividual);
subplot(Ndata,1,k);
h1 = 0;
if strcmp(data(1).Type,'Cone')
h1=plot(data(k).T, data(k).M(:,2),'k');  
else
h1=plot(data(k).T, data(k).M,'k');
end
hold on;
h2=plot(T(:,2),M(:,2), 'k');

set(h2, 'LineStyle', '--');

if strcmp(data(1).Type,'Cone') && strcmp(data(1).gas, 'Air')
set(0,'CurrentFigure', bestIndividual_2);
subplot(Ndata,1,k);
h1=plot(data(k).T, data(k).M(:,1),'k');
hold on;
h2=plot(T(:,2),M(:,1), 'k'); 
set(h2, 'LineStyle', '--');
end

end %k

if strcmp(data(1).Type,'TGA')
xlabel('Temperature [\circC]');
ylabel('Mass fraction');
elseif strcmp(data(1).Type,'Cone')
    set(0,'CurrentFigure', bestIndividual);
    xlabel('Time [s]');
    ylabel('MLR [kg/m^2s]');
    set(0,'CurrentFigure', bestIndividual_2);
    xlabel('Time [s]');
    ylabel('HRR [kW/m^2]');
    set(bestIndividual_2, 'Visible', 'on');
end
set(bestIndividual, 'Visible', 'on');

gen = gen+1;
end %end of while

%---------------------- GA loop ends here

%plot best individual
[oV, line] = min(objVal);
bestChrom = Chrom(line,:);
if ishandle(bestIndividual)
    delete(bestIndividual);
    bestIndividual = figure('Visible', 'off', ...
                    'Name', 'Best Individual', ...
                    'NumberTitle', 'off');
end

if ishandle(bestIndividual_2)
    delete(bestIndividual_2);
    bestIndividual_2 = figure('Visible', 'off', ...
                    'Name', 'Best Individual', ...
                    'NumberTitle', 'off');
end

for k=1:Ndata
    
    if strcmp(data(1).Type,'TGA')
        parameters = zeros(1,length(parInd));
        index = 5;
        NR = par(5);
        L = length(par);
        for n=1:length(parInd)
    
         if parInd(n)==1
             %if tramp T parameter
                if par(1)==1
                parameters(n) = par(2)/data(k).Rate*60;
                %if residue parameter
                %if mass fraction parameter
                elseif parInd(n)==length(parInd)
                 parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L))); 
                 index = index+1;
                 %if residue parameter
                else
                index = index+1;
                parameters(n)=1-Chrom(i,par(index));
                end
        elseif parInd(n)==2
                if par(1)==1 && par(3)==1
                 parameters(n)= par(4)/data(k).Rate*60;
                elseif parInd(n)==length(parInd)
                 parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L))); 
                 index = index+1;
                 %if residue parameter
                else
                index = index+1;
                parameters(n)=1-Chrom(i,par(index));
                end
        elseif parInd(n)==length(parInd) && parInd(n) > (2+par(5))
                parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L)));
                index = index +1;
        else
            index = index+1;
            parameters(n)=1-Chrom(i,par(index));
            %if mass fraction parameter
     
        end%end of if parInd(n)
        end %end of for n
    
    [M,T]=fds_tga(template,FdsExe,'tga',bestChrom,LogScaling,parameters,[],0);
    M(:,2)=M(:,2)./max(M(:,2));
    dt = T(2,1)-T(1,1); %time step
    T0 = T(1,2);
       if max(T(:,2))<max(data(k).T)
            i=length(T(:,2));
            while T(i-1,2)<max(data(k).T)
                M(i,1)=M(i-1,1)+dt;
                M(i,2)=M(i-1,2);
                T(i,1)=T(i-1,1)+dt;
                T(i,2)=(T(i,1)-dt)*data(k).Rate/60+T0;
                i=i+1;
            end
        end
   
elseif strcmp(data(1).Type,'Cone')
    if strcmp(data(k).gas,'Air')
        temp = template;
    else %if in N2
        temp = 'ox_limited_input.fds';
    end
parameters = [];
[M,T]=fds_tga(temp,FdsExe,'cone',bestChrom,LogScaling,parameters,[],data(k).dataType);
M(:,2)=coneFilter(T(:,2),M(:,2));
M(:,1)=M(:,1)./data(1).A;
end

set(0,'CurrentFigure', bestIndividual);
subplot(Ndata,1,k);
h1 = 0;
if strcmp(data(1).Type,'Cone')
h1=plot(data(k).T, data(k).M(:,2),'k');  
else
h1=plot(data(k).T, data(k).M,'k');
end
hold on;
h2=plot(T(:,2), M(:,2), 'k');
 
set(h2, 'LineStyle', '--');

if strcmp(data(1).Type,'Cone')&& strcmp(data(1).gas, 'Air')
    set(0,'CurrentFigure', bestIndividual_2);
    subplot(Ndata,1,k);
    h1=plot(data(k).T, data(k).M(:,1),'k');
    hold on;
    h2=plot(T(:,2),M(:,1), 'k'); 
    set(h2, 'LineStyle', '--');
    set(0,'CurrentFigure', bestIndividual_2);
    xlabel('Time [s]');
    ylabel('HRR [kW/m^2]');
    set(bestIndividual_2, 'Visible', 'on');
end
end % end of for k

if strcmp(data(1).Type,'TGA')
    xlabel('Temperature [\circC]');
    ylabel('Mass fraction');
elseif strcmp(data(1).Type,'Cone')
    set(0,'CurrentFigure', bestIndividual);
    xlabel('Time [s]');
    ylabel('MLR [kg/m^2s]');   
end

if strcmp(data(1).Type,'TGA')
%plot gradient
parameters = zeros(1,length(parInd));
index = 5;
NR = par(5);
    L = length(par);
    for n=1:length(parInd)
      if parInd(n)==1
         %if tramp T parameter
            if par(1)==1
            parameters(n) = par(2)/data(k).Rate*60;
            %if residue parameter
            %if mass fraction parameter
            elseif parInd(n)==length(parInd)
             parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L))); 
             index = index+1;
             %if residue parameter
            else
            index = index+1;
            parameters(n)=1-Chrom(i,par(index));
            end
        elseif parInd(n)==2
            if par(1)==1 && par(3)==1
             parameters(n)= par(4)/data(k).Rate*60;
            elseif parInd(n)==length(parInd)
             parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L))); 
             index = index+1;
             %if residue parameter
            else
            index = index+1;
            parameters(n)=1-Chrom(i,par(index));
            end
        elseif parInd(n)==length(parInd) && parInd(n) > (2+par(5))
            parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L)));
            index = index +1;
        else
            index = index+1;
            parameters(n)=1-Chrom(i,par(index));
            %if mass fraction parameter

        end
    end %end of for n
[M,T]=fds_tga(template,FdsExe,'tga',bestChrom,LogScaling,parameters,[],0);
M(:,2)=M(:,2)./max(M(:,2));
d = M(2,1)-M(1,1);
g = -gradient(M(:,2),d);
t=M(:,1);
grad = figure('Visible', 'off', ...
                    'Name', 'Gradient', ...
                    'NumberTitle', 'off');
plot(data(1).Time, -data(1).dMdt, t, g);
legend('Exp', 'Model');
xlabel('Time [s]');
ylabel('dM/dt [1/s]');
set(grad, 'Visible', 'on');
saveas(grad, 'gradient.fig', 'fig');
end % end of if TGA

set(bestIndividual, 'Visible', 'on');

%save best input file
filename = create_input(template, bestChrom, LogScaling,parameters);
msgbox(['Input file saved as ' filename]);
save('best.mat', 'bestChrom');
saveas(bestIndividual,'bestInd.fig','fig');
saveas(fitness,'fit.fig','fig');
if exist('ga.stop','file')
    delete ga.stop
end
bestChrom(iLogScaling) = 10.^bestChrom(iLogScaling);
end %end of function

%----------------------------------------------------------------------

function objVal=objF(Ndata, data, Chrom, LogScaling, template, FdsExe, weights, Aindex, limits, par, parInd, moisture)
%time, data are experimental data
    
    %weights = [985/1000 15/1000]; %weights of fitness function: 1. data 2. gradient 3.penalty of A
    [NIND,NVAR] = size(Chrom);  
    R=[];
    R_1 = [];
    R_2 = [];
    for i=1:NIND
    SV = 0; %SV is the scaling value, sum of squares
    rdata = [];
    dataMod = [];
    for k = 1:Ndata
       if strcmp(data(1).Type,'TGA')
        parameters = zeros(1,length(parInd));
        index = 5;
        NR = par(5);
        L = length(par);
       for n=1:length(parInd) 
           
            if parInd(n)==1
                 %if tramp T parameter
                    if par(1)==1
                    parameters(n) = par(2)/data(k).Rate*60;
                    %if residue parameter
                    %if mass fraction parameter
                    elseif parInd(n)==length(parInd)
                     parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L))); 
                     index = index+1;
                     %if residue parameter
                    else
                    index = index+1;
                    parameters(n)=1-Chrom(i,par(index));
                    end
            elseif parInd(n)==2
                    if par(1)==1 && par(3)==1
                     parameters(n)= par(4)/data(k).Rate*60;
                    elseif parInd(n)==length(parInd)
                     parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L))); 
                     index = index+1;
                     %if residue parameter
                    else
                    index = index+1;
                    parameters(n)=1-Chrom(i,par(index));
                    end
            elseif parInd(n)==length(parInd) && parInd(n) > (2+par(5))
                    parameters(n)=1-moisture-sum(Chrom(i,par(6+NR:L)));
                    index = index +1;
            else
                index = index+1;
                parameters(n)=1-Chrom(i,par(index));
                %if mass fraction parameter
     
            end
       
       end % end of for n
       end %if TGA
     resMass = min(data(k).M);
        if strcmp(data(1).Type,'TGA')
        [M,T]=fds_tga(template,FdsExe,'tga',Chrom(i,:),LogScaling,parameters,[],0);
        elseif strcmp(data(1).Type,'Cone')
        if strcmp(data(k).gas,'Air')
            temp = template;
        else %if in N2
            temp = 'ox_limited_input.fds';
        end
        parameters = [];
        [M,T]=fds_tga(temp,FdsExe,'cone',Chrom(i,:),LogScaling,parameters,[],data(k).dataType);
        M(:,1)=coneFilter(T(:,2),M(:,1));
        M(:,2)=coneFilter(T(:,2),M(:,2));
        M(:,1)=M(:,1)./data(1).A;
        end
     fclose('all');
     l=1;
     if min(T(:,2))>min(data(k).T)        
            while min(T(:,2))>data(k).T(l)
                l=l+1;
            end
     end
    Mass = [];
    Mass(:,1) = interp1(T(:,2), M(:,1), data(k).T(l:length(data(k).T)));
    Mass(:,2) = interp1(T(:,2), M(:,2), data(k).T(l:length(data(k).T)));
    if isnan(Mass(1))
        display('nyt');
    end
    
    Mass(:,1) = removeNaNM(Mass(:,1),1);
    Mass(:,2) = removeNaNM(Mass(:,2),1);
    
    if strcmp(data(1).Type,'TGA')
    %Mass
    dataMod = data(k).M(l:length(data(k).T));   
    Mass(:,2)=Mass(:,2)./max(Mass(:,2)); %scale
    SV = sum((dataMod-mean(dataMod)).^2);   
    R(k) = sum((dataMod-Mass(:,2)).^2); 
    rdata_mass = weights(1)*(1-(SV-R(k))/SV);
    
    %Gradient
    grad_Mod = -gradient(data(k).M(l:length(data(k).T)),data(k).T(l:length(data(k).T)));   
    grad_exp = -gradient(Mass(:,2),data(k).T(l:length(data(k).T)));
    SV_grad = sum((grad_Mod-mean(grad_Mod)).^2);   
    R_grad = sum((grad_Mod-grad_exp).^2); 
    rdata_grad = weights(2)*(1-(SV_grad-R_grad)/SV_grad);
    
    rdata(k) = rdata_mass + rdata_grad;
    
    elseif strcmp(data(1).Type,'Cone')     
    dataMod(:,1) = data(k).M(l:length(data(k).T),2);     
    SV_2 = sum((dataMod(:,1)-mean(dataMod(:,1))).^2);  
    R_2(k) = sum((dataMod(:,1)-Mass(:,2)).^2); 
    rdata_2 = 1-(SV_2-R_2(k))/SV_2;
    
    if strcmp(data(k).gas,'Air')
            dataMod(:,2) = data(k).M(l:length(data(k).T),1);
            SV_1 = sum((dataMod(:,2)-mean(dataMod(:,2))).^2); 
            R_1(k) = sum((dataMod(:,2)-Mass(:,1)).^2); 
            rdata_1 = 1-(SV_1-R_1(k))/SV_1;
            rdata(k) = (rdata_1+ rdata_1)/2;
            
    else
        rdata(k) = rdata_2;
        
    end
    end
    
    end % end of for k
     A=0;
    if any(Aindex)  
    for j = 1:length(Aindex)
        A = A + (Chrom(i,Aindex(j))-limits(1,Aindex(j)))./(limits(2,Aindex(j))-limits(1, Aindex(j)));
    end
    A = A/length(Aindex(j));
    end % end of if any
    % final objective value is weighted mean of models with different rates
    objVal(i) =mean(rdata)+weights(3)*A;
    end %end of for i
    objVal = objVal';
    
end % end of function


