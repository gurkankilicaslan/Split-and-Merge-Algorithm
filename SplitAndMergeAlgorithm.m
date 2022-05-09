clear all; close all; clc;


%% PART I
load("Data.mat")


mean_range = mean(rangedata_Q1); %Take average of range data

%Calculate x and y values:
for i=1:length(rangedata_Q1)
    x(i) = rangedata_Q1(i)*cosd(angledata_Q1(i));
    y(i) = rangedata_Q1(i)*sind(angledata_Q1(i));
    weight(i) = 1/(rangedata_Q1(i))^2; %weight calculation
end

sum_weight = sum(weight); %sum of weights

% Calculating alpha formula which includes single and double summations
%elements:
a = 0;
for i=1:length(rangedata_Q1)
    a = a + (weight(i)*rangedata_Q1(i)^2)*sind(2*angledata_Q1(i));
end

b = 0;
for i=1:length(rangedata_Q1)
    for j=1:length(rangedata_Q1)
        b = b + weight(i)*weight(j)*rangedata_Q1(i)*rangedata_Q1(j)...
            *cosd(angledata_Q1(i))*sind(angledata_Q1(j));
    end
end

c = 0;
for i=1:length(rangedata_Q1)
    c = c + (weight(i)*rangedata_Q1(i)^2)*cosd(2*angledata_Q1(i));
end

d = 0;
for i=1:length(rangedata_Q1)
    for j=1:length(rangedata_Q1)
        d = d + weight(i)*weight(i)*rangedata_Q1(i)*rangedata_Q1(j)...
            *cosd(angledata_Q1(i) + angledata_Q1(j));
    end
end

% alpha calculation at the end:
alpha_rad = 0.5*atan((a-(2/sum_weight)*b)/(c-(1/sum_weight)*d));

alpha = 180*alpha_rad/pi; %radian to degrees

% r calculation:
r1 = 0;
for i=1:length(rangedata_Q1)
    r1 = r1 + weight(i)*rangedata_Q1(i)*cosd(angledata_Q1(i)-alpha);
end
r = r1/sum_weight;

fh = figure();
fh.WindowState = 'maximized';

scatter(x,y,20,"k")
hold on
hold on


x1 = r/cosd(alpha); %intersection of fitted line and x axis
y1 = r/sind(alpha); %intersection of fitted line and y axis

% Fitted Line parameters:
x_line = -3:0.1:1.5;
y_line = tand(90+alpha)*x_line + y1;

plot(x_line,y_line,'LineWidth',2,color = "r")
title("Fitted Line of Sensor Data")
xlabel("X")
ylabel("Y")
hold on
text(0.58,-0.085,"\leftarrow fitted line",'FontSize',12)

%Standard error of the mean

x_error = std(x);
y_error = std(y);
for i=1:length(rangedata_Q1)
    standard_error_for_x(i) = x_error/sqrt(length(x));
    standard_error_for_y(i) = y_error/sqrt(length(y));
 end

e = errorbar(x,y,standard_error_for_x,'o'); %error bar for x
e.Color = 'black';

e = errorbar(x,y,standard_error_for_y,'horizontal','o'); %error bar for y
e.Color = 'black';



% End points of the r vector
x_c = r*cosd(alpha);
y_c = r*sind(alpha);

% Draw r vector
quiver(0,0,x_c,y_c,0,'LineWidth',2,'Color','k')
text(x_c/2,y_c/2,"\leftarrow r = " + string(r),'FontSize',12)

axis([-0.05 0.6 -0.15 0.65])
axis equal

% Draw angle alpha
th = 0:1/360:alpha*pi/180;  radi = 0.1;
xx = radi*cos(th); yy = radi*sin(th);
plot(xx,yy,'LineWidth',2,'Color','b'),grid
text(0.1,0.03,'\leftarrow \alpha = '+string(alpha) + " deg.",'FontSize',12)

legend("sensor data","fitted line","y_{error bars}", "x_{error bars}",...
    "r vector")



%% PART II

%Functions created for this part(Question 2) as follow:
%Point_Line_Distance --> it finds distance between point and line
%Difference_2points --> it finds distance between 2 points
%point_to_line --> it finds distance between point and line, different use
% You can find them at the end of the code


Threshold_For_Error = 0.2; %threshold kinda says where to stop splitting
% You can change the threshold

load("Data.mat")
for i=1:length(rangedata_Q2)
    x(i) = rangedata_Q2(i)*cosd(angledata_Q2(i));
    y(i) = rangedata_Q2(i)*sind(angledata_Q2(i));
end
DataMatrix(:,1) = x; %Filling general Data Matrix
DataMatrix(:,2) = y; %Filling general Data Matrix

fit_prop = []; %Fitting matrix, to be filled

DistanceMax = 0; %To keep track of where to split lines

%Initialize Line and Distance Values
for i=1:1:size(DataMatrix,1)
    DataMatrix(i,3) = 1; %What line that point belongs to
    DataMatrix(i,4) = 1000; % 1000 is arbitrarily, just to be high
    DataMatrix(i,5) = 0; % it is 0 if it is an end point, 1 if not
end

%Find which points are the furthest away
for i = 1:1:size(DataMatrix,1)       
    for j = 1:1:size(DataMatrix,1)
        Distance_Now = Difference_2points(DataMatrix(i,1), ...
            DataMatrix(j,1),DataMatrix(j,2)); %calculates current distance
        if (Distance_Now > DistanceMax) %if this occurs, it keeps iterating
            DistanceMax = Distance_Now;
            fit_prop(1) = DataMatrix(i,1);
            fit_prop(2) = DataMatrix(i,2);
            DataMatrix(i,5) = 1; %It is an end point
            fit_prop(3) = DataMatrix(j,1);
            fit_prop(4) = DataMatrix(j,2);
            DataMatrix(j,5) = 1; %It is an end point
        end
    end
end
%Initialize distance from each point to the line

Error_Huge = 1000; % 1000 is just to be huge, to enter the while loop

%Splitting Segments starts here
while (Error_Huge >= Threshold_For_Error)
    
    %Check to see if any of the distances exceed the Distance.
    Error_Huge = 0; % Since we are in, it is 0 now
    Out = 1; %Outlier points
    for i = 1:1:size(DataMatrix,1)
      DataMatrix(i,4) = Point_Line_Distance(fit_prop(DataMatrix(i,3),:),...
            DataMatrix(i,1),DataMatrix(i,2)); %Find distance point to line
        
        if  (DataMatrix(i,4) > Error_Huge)
            %If we've arrived at this point, we should see 
            %if this point belongs to a different line.
            for j = 1:1:size(fit_prop,1)
              Error = Point_Line_Distance(fit_prop(j,:),DataMatrix(i,1),...
                    DataMatrix(i,2)); %Find distance point to line
                
              if (Error < DataMatrix(i,4))
                    DataMatrix(i,3) = j; %Change the line point belongs to
                    DataMatrix(i,4) = Error; %Change the error
              end
            end
            %If we still see a big error;
            if (DataMatrix(i,4) > Error_Huge)
                Error_Huge = DataMatrix(i,4);
                Out = i; %Outlier Points
            end
        end
    end
    
    
    %It is necessary to break the line in which the outlier is a member.
    if (Error_Huge >= Threshold_For_Error)
        Divide_Line = DataMatrix(Out,3); %Line we now split.
        DataMatrix(Out,5) = 1; %It is an end point
        X_Created = DataMatrix(Out,1); %New x value
        Y_Created = DataMatrix(Out,2); %New y value
        
        %New line segment:
        %Update fitting matrix        
        fit_prop = [fit_prop; X_Created Y_Created ...
            fit_prop(Divide_Line, 3) fit_prop(Divide_Line,4)]; 
        
        %Update the original segment
        fit_prop(Divide_Line,3) = X_Created;
        fit_prop(Divide_Line,4) = Y_Created;
    end
    
end
fh = figure();
fh.WindowState = 'maximized';

plot(DataMatrix(:,1),DataMatrix(:,2),'.','MarkerSize',10) %Plot sensor data
title({'Split and Merge Method for given Sensor Data',...
    'Threshold you chose = ' + string(Threshold_For_Error)}) 
xlabel("X")
ylabel("Y")
axis([-4 4 -3.2 4.2])
grid on
hold on

for i = 1:1:size(fit_prop,1)
    %Plot fitted lines:
    plot([fit_prop(i,1) fit_prop(i,3)], [fit_prop(i,2) fit_prop(i,4)],...
        'Color','r','LineWidth',1.2);
    
    %Calculate r(distance from (0,0) to fitted line):
    r = point_to_line([0,0,0], [fit_prop(i,1),fit_prop(i,2),0], ...
        [fit_prop(i,3),fit_prop(i,4),0]);
    
    %Calculate angle between x axis and r vector:
    alpha = acotd((fit_prop(i,2) - fit_prop(i,4))/(fit_prop(i,3) ...
        - fit_prop(i,1)));
    
    %Some necessary changes to draw r vector:
    if fit_prop(i,2)<0 && fit_prop(i,3)<0 && fit_prop(i,4)<0 ...
            && fit_prop(i,1)<0
        r=-r;
    end
    
    if fit_prop(i,1) <0 || fit_prop(i,3)<0
        r=-r;
    end
    
    if fit_prop(i,2) <0 && fit_prop(i,4)<0
        r=-r;
    end
    
    %Calculate end points of r vector:
    x_c = r*cosd(alpha);
    y_c = r*sind(alpha);
    
    %Draw r vector
    quiver(0,0,x_c,y_c,0,'LineWidth',1,'Color','k')
    text(x_c/2, y_c/2, "r = "+r + ", \alpha = " + alpha,'FontSize',8)
    
end
axis equal
legend("sensor data","fitted lines","r vectors")


%Following functions are used throughout the code. What they do is 
%mentioned at the very beginning of Question 2 part:

function Distance = Point_Line_Distance(line_to_point,x_0,y_0)
    x_1 = line_to_point(1);
    y_1 = line_to_point(2);
    x_2 = line_to_point(3);
    y_2 = line_to_point(4);
    
    %Determine the two orthogonal lines
    %Find the slope of the line to point
    slope = (y_2 - y_1)/(x_2 - x_1);
    %Now we need the inverse reciperocal
    slope = -1/slope;
    if (slope == -inf)||(slope == inf)
        slope = 100000; %100000 is arbitrary, just to be huge
    end
    a_1 = slope;
    b_1 = -1;
    c_1 = y_1 - slope*x_1;
    a_2 = slope;
    b_2 = -1;
    c_2 = y_2 - slope*x_2;
    
    d_1 = abs(a_1*x_0+b_1*y_0+c_1)/sqrt(a_1^2+b_1^2);
    d_2 = abs(a_2*x_0+b_2*y_0+c_2)/sqrt(a_2^2+b_2^2);
    d_0 = sqrt((x_2 - x_1)^2 + (y_2 - y_1)^2);
    if ((d_1 > d_0) || (d_2 > d_0))
        Distance = min(sqrt((x_2-x_0)^2+(y_2-y_0)^2), ...
            sqrt((x_0-x_1)^2+(y_0-y_1)^2));
    else
        numerator = abs((y_2-y_1)*x_0-(x_2-x_1)*y_0+x_2*y_1-y_2*x_1);
        denominator = sqrt((y_2-y_1)^2+(x_2-x_1)^2);
        Distance = numerator/denominator;
    end
end

function distance = Difference_2points(a,b,c) % distance finder
    distance = sqrt((b-c)^2 + (a-c)^2);    
end

function d = point_to_line(pt, v1, v2) % r finder
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end
