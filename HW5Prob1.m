
tmp        = csvread('rates.csv',1); %ignore first row, since it has text column labels
GPS_lon    = tmp(:,1); %Longitude
GPS_lat    = tmp(:,2); %Latitude
GPS_x      = tmp(:,3); %Easting (in meters, UTM zone 11)
GPS_y      = tmp(:,4); %Northing (in meters, UTM zone 11)
GPS_UE     = tmp(:,5); %Displacement rate in the east direction, mm/yr 
GPS_UN     = tmp(:,6); %Displacement rate in the north direction, mm/yr

%faults and coastlines for plotting.
tmp        = csvread('faults.csv',1);
faults_lon = tmp(:,1);
faults_lat = tmp(:,2);
faults_x =   tmp(:,3); 
faults_y =   tmp(:,4);

tmp        = csvread('coasts.csv',1);
coast_lon  = tmp(:,1);
coast_lat  = tmp(:,2);
coast_x    = tmp(:,3); %added UTM coordinates to the coasts for plotting
coast_y    = tmp(:,4);

%the original "map" locations have values like 4e6 meters - this just makes
%the plots a little cleaner.  Doesn't affect any of the math at all.
coast_x   = coast_x-mean(GPS_x);
coast_y   = coast_y-mean(GPS_y);
faults_x  = faults_x-mean(GPS_x);
faults_y  = faults_y-mean(GPS_y);
GPS_x     = GPS_x-mean(GPS_x);
GPS_y     = GPS_y-mean(GPS_y);

GPS_UE     = GPS_UE/1000;  %convert to meters/yr, since our units of distance are in meters
GPS_UN     = GPS_UN/1000;

%%%%
%%%% All the stuff above is old, the stuff below is new. Note how figures
%%%% have "if" statements so that they only show up if you change the 0 to
%%%% a 1. This is to help you walk through the code and make changes. Note
%%%% that you should "close all" from the command line if you have too many
%%%% figures open.

%%% Plot velocity vectors on normal plot (E,W axes)
%No changes necessary here, but read notes.
if(1)
    %Note that when we remove the mean of the GPS vectores, it doesn't
    %necessarily mean that all of the vectors on one side of the fault go
    %one way and the ones on the other side go the other way.  If we have
    %more points on one side, then you'll get a mean value weighted towards
    %those values.  But it makes it a little easier to see the complexity.
    figure
    plot(coast_x/1e3,coast_y/1e3,'b')
    hold on
    plot(faults_x/1e3,faults_y/1e3,'r')
    quiver(GPS_x/1e3,GPS_y/1e3,GPS_UE-mean(GPS_UE),GPS_UN-mean(GPS_UN),'k')
    axis image
    axis([-3.5e2 3e2 -2.5e2 2e2])
    xlabel('Easting (km)')
    ylabel('Northing (km)')
end    


%%%Rotate region
%We want to plot the rates along a profile perpendicular to the fault.  One
%way to do this is to rotate the whole region, and then just pull out a
%given range of "Y" values in the new coordinate system.  Try various
%values of "theta", for the angle of rotation, until the San Andreas is
%essentially up-down on your plot (i.e., vectors on one side point up, on
%the other side they point down).
%%% CHANGE THIS
theta = 10; %in degrees

%This is a 2x2 matrix to rotate the coordinates of coastlines, faults, etc.
rotM = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];

%Now rotate all the things we've been plotting
M     = rotM*[faults_x faults_y]';
rotfx = M(1,:);
rotfy = M(2,:);

M   = rotM*[coast_x coast_y]';
rotcx = M(1,:);
rotcy = M(2,:);

M =  rotM*[GPS_x GPS_y]';
rotGPS_x = M(1,:);
rotGPS_y = M(2,:);

M =  rotM*[GPS_UE GPS_UN]';
rotGPS_UE = M(1,:);
rotGPS_UN = M(2,:);

if(0)
    %Plot the variables in the rotated coordinate system
    figure
    quiver(rotGPS_x,rotGPS_y,rotGPS_UE-mean(rotGPS_UE),rotGPS_UN-mean(rotGPS_UN),'k')
    axis image
    a=axis;
    hold on
    plot(rotfx,rotfy,'r')
    plot(rotcx,rotcy,'b')
    axis(a)
end


%%% Select range of Y-values and plot profile
% The San Andreas does not run perfectly straigh - it has a stepover north
% of Los Angeles.  In your rotated plot, you should be able to see that
% there is a horizontal offset between the upper and lower sections.  If
% you plot the whole thing vs your new, rotated X coordinate, you will get
% the following plot.

if(0)
    figure
    plot(rotGPS_x,rotGPS_UN,'r.')
    xlabel('Distance across fault, rotated coords, meters')
    ylabel('Fault-parallel displacement rate (m/yr)')
    grid on
end


% the deformation zone looks really broad and smeared out because you are
% combining the northern and southern sections together.  Instead, let's
% pick JUST the southern section.  To do that, we will use ONLY values that
% fall within certain bounds in the rotated Y variable. Use the Y axis as
% shown in your previous, rotated map view plot.

%%%CHANGE THESE to actual numbers, they are set to a default max/min in your rotated coords.
min_roty = min(rotGPS_y);
max_roty = max(rotGPS_y);
  
goodid = find(rotGPS_y < max_roty & rotGPS_y > min_roty);

if(0)
    figure
    hold on
    plot(rotGPS_x(goodid)/1e3,rotGPS_UN(goodid),'r.')
    xlabel('Distance across fault, rotated coords, km')
    ylabel('Fault-parallel displacement rate (m/yr)')
    grid on
end


%%% Now you need to model the fault, using the equation shown from class. Tinker with the following

%%% Location of fault (this should be about where the center of your
%%% "arctan" signal is, in meters (note that the last figure shows km)
x0 = 5e3;

%%% Locking depth, in meters
D = 35e3;

%%% Plate rate (in meters/ year, like your UE, UN vectors.  Should be about
%%% the total from one side of the fault to the other. The negative sign
%%% makes the value right-lateral (left side positive)
V = -0.08;

%Predicted profile of displacement rate, at just our selected points.  
%%%This equation is WRONG!  Fix it to what we showed in class.
model = V/pi*acos(D/(rotGPS_x(goodid)-x0));

%solve for the average offset of this curve (i.e., center on the data)
meanrate = mean(rotGPS_UN(goodid)-model);
  
%std.dev of residual
rms = std(rotGPS_UN(goodid)-model);

if(0)
    %plot results!
    figure
    plot((rotGPS_x(goodid)-x0)/1e3,rotGPS_UN(goodid)-meanrate,'.')
    hold on
    plot((rotGPS_x(goodid)-x0)/1e3,model,'r.')
    grid on
    xlabel('Distance across fault, rotated coords, km')
    ylabel('Fault-parallel displacement rate (m/yr)')
    legend('Data','Model')
    title(['std dev of res. = ' num2str(rms*100) ' cm/yr'])
end
  

