%% ELEC 4700 Assignment 3 Question One Work

% Name: Seth Thompson
% Student Number: 101031310

close all
clear
clc

% 1 A) A voltage of 0.1V is applied across the x-dimension, what is the
% electric field on the electrons? Voltage is constant over the region.

% Leftside of the region is at 0.1V, right side is at 0V, Linspace can be
% used to make a matrix with the appropriate values

regionLength = 200e-9; % meters
regionWidth = 100e-9; % meters

% Create a vector to act as the voltage drop from 0.1V to 0V over the
% region

voltageVector = linspace(0.1,0,200);

% Copy the elements of the voltage vector into each row of a voltage matrix

for n = 1:100
    voltageMatrix(n,:) = voltageVector;
end

% Plotting the voltage over the region

figure(1)
surf(voltageMatrix)
title({'Voltage Over the Semiconductor','Seth Thompson | 101031310'})
xlabel('Semiconductor Length')
ylabel('SemiConductor Width')
zlabel('Voltage (V)')
grid on

pause(1)

% The electric field over the region can be obtained by taking the gradient
% of the negative voltage drop (E = -del(V))

[ElecFieldX,ElecFieldY] = gradient(-voltageMatrix,1e-9);

% Plotting the Electric field over the region

figure(2)
quiver(ElecFieldX,ElecFieldY)
title({'Electric Field Over the Semiconductor','Seth Thompson | 101031310'})
xlabel('Semiconductor Length')
ylabel('SemiConductor Width')
xlim([0,regionLength*1e9])
ylim([0,regionWidth*1e9])

pause(1)

% To get the value of the electric field at a certain point, use the
% equation del(V) = -Ed. In this case del(V) is 0.1V and d is 200nm
% Outputting the constant value of the electric field.

EFieldVal = ElecFieldX(1,1);
EFieldValDISPLAY = EFieldVal*1e-3;
fprintf('The electric field is constant over the semiconductor.\n')
fprintf('It has no component in the y-direction and has a value of %f kV/m in the x direction.\n\n',EFieldValDISPLAY)
% These results make sense, the voltage is only changing along the
% x-direction and not the y-direction, so taking the gradient of the
% voltage would result in an electric field going from left to right along
% the x-axis and no electric field along the y-axis.

% On top of that, since the voltage drop from the left to the right drops
% linearlly, the resulting electric field should be a constant. If the
% gradient takes the slope over the region, and the slope at certain points
% of the region is always constant, then it would make sense to have the
% result be a constant.

% 1 B) What is the force on each electron?

% It is known that F = Eq, so multiplying the electric field by the
% effective charge will give the force on each electron.

eCharge = 1.602e-19; % Columbs
Force = EFieldVal*eCharge; % Newtons
ForceDISPLAY = Force*(10^15);
fprintf('With respect to the force on each electron, there is no component\n')
fprintf('in the y-direction, but in the x-direction it has a value of %f fN.\n\n',ForceDISPLAY)

% 1C) Calculate the acceleration on the electrons and use this in your model to
% update the velocity of each electron at each time step. Add the capability
% for the electrons to respond to a static electric field with both an x 
% and a y component. Plot the trajectories of the electrons. They should be
% curved! Increase the electric field and see the curve! Be careful that 
% your time step is appropriate!

% To calculate the acceleration on the electrons, the equation F = m*a can
% be rearranged to solve for the acceleration since we know both the force
% and the mass. This is shown next

mRest = 9.109e-31; % kilograms
mEffective = 0.26*mRest; % kilograms
acceleration = Force/mEffective;
accelerationDISPLAY = acceleration*1e-15;
fprintf('With respect to the acceleration on each electron, there is no component\n')
fprintf('in the y-direction, but in the x-direction it has a value of %f Pm/s^2.\n\n',accelerationDISPLAY)

% Code from assignment 1 will be used again to make the 2D plot of particle
% trajectories.

% Defining Constants to be used in this part of the assignment.
regionLength = 200e-9; % meters
regionWidth = 100e-9; % meters
Temperature = 300; % Kelvin
kb = 1.380649e-23; % J*K^-1

% Equation for thermal velocity is sqrt((2*kb*T)/m)(RMS)

thermalVel = sqrt((2*kb*Temperature)/mEffective);

% Defining the given constant
meanTime = 0.2e-12; % seconds

% Assigning each particle a random starting position in the 100x200 plane.

electronAmount = 10; % Amount of electrons being simulated

% Assigning a particle a random position on the XY plane within the maximum
% and minimum width and length.
particleXPosition = regionLength.*rand(electronAmount,2);
particleYPosition = regionWidth.*rand(electronAmount,2);

% First column is the initial position, second column is the next position,
% think xn and xn+1
particleXPosition(:,1) = particleXPosition(:,2);
particleYPosition(:,1) = particleYPosition(:,2);

% Timestep calculated from sqrt(100nm^2 + 200nm^2)/1000, which is approx
% 0.22nm, and then 0.22nm/thermalVel
timeStep = 0.22e-9/thermalVel; % seconds

% Assigning a particle a random velocity using the randn function
% (normally distributed random values).
% ** THIS TIME, the velocity from the acceleration in the x-direction must
% be added **
randXV = randn(electronAmount,2)*sqrt((kb*Temperature)/mEffective);
randYV = randn(electronAmount,2)*sqrt((kb*Temperature)/mEffective);

% Creating a loop that will calculate each particles position as time goes
% on. Also making a counter that will keep track of the amount of
% iterations done. 
counter = 0;

% Calculating the displacement of the electron for each time step
particleXDisplacement = timeStep*randXV(:,1);
particleYDisplacement = timeStep*randYV(:,1);

% Creating a vector of random RGB values so each particle will have its own
% colour
colourVector = rand(electronAmount,3);

% Creating avariable to hold the probabillity of an electron scattering
PScat = 1 - exp(-(timeStep/meanTime));

% Creating two vectors that will hold the initial position before a particle
% has deflected. One for the X position and another for the Y. There will
% also be a collide counter to determine the aount of collisions made
% during the system

beforeCollideX = particleXPosition(:,1);
beforeCollideY = particleYPosition(:,1);

for m = 1:500
    % This time, velocities must be updated over time because of the
    % acceleration now present in the x-direction.
    randXV = randXV + acceleration*timeStep;
    % Updating the position of the particle
    for n = 1:electronAmount
        % First checking to see if the electron scatters, if it does then a
        % new velocity must be generated and from there a new displacement
        if (rand < PScat)
            thetaScatter = (2*pi).*rand(1);
            randXV(n,:) = randn(1)*sqrt((kb*Temperature)/mEffective) + acceleration*timeStep;
            randYV(n,:) = randn(1)*sqrt((kb*Temperature)/mEffective);
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleYDisplacement(n,:) = timeStep*randYV(n,1);           
            
            % Saving the new positions before the next collision
            beforeCollideX = particleXPosition(:,1);
            beforeCollideY = particleYPosition(:,1);
        end
        % Checking the boundary conditions for the left and right sides of
        % the plot. If the particle is going to pass a border on the left
        % or right side then simply move it to the left or right.
        if (particleXPosition(n,1) + particleXDisplacement(n) > 2e-7)
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1) - 2e-7;
        elseif (particleXPosition(n,1) + particleXDisplacement(n) < 0)
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1) + 2e-7;
        else
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1);
        end
        
        % Checking boundary conditions for the top and bottom of the plot.
        % If the particle is about to cross the top of the plot, flip the
        % y, displacement value and have it move in the opposite direction;
        % the opposite applies for the bottom of the plot.
        if ((particleYPosition(n,1) + particleYDisplacement(n) > 1e-7) || (particleYPosition(n,1) + particleYDisplacement(n) < 0))
            particleYDisplacement(n) = -particleYDisplacement(n);
            particleYPosition(n,2) = particleYDisplacement(n) + particleYPosition(n,1);
        else
            particleYPosition(n,2) = particleYDisplacement(n) + particleYPosition(n,1);
        end
    end
    
    % Creating a vector that will hold all timesteps over the simulation
    timeVector(m) = timeStep*m;
    
    if (counter == 0)
        % Plotting the particles over time.
        figure(3)
        scatter(particleXPosition(:,2),particleYPosition(:,2),1,colourVector(:,1))
        hold on
        title({['2-D Plot of Particle Trajectories'],['Seth Thompson | 101031310']})
        xlabel('X-Axis (m)')
        ylabel('Y-Axis (m)')
        xlim([0 200e-9])
        ylim([0 100e-9])
    elseif (counter > 0 && counter < 1000)
        title({['2-D Plot of Particle Trajectories'],['Seth Thompson | 101031310']})
        scatter(particleXPosition(:,2),particleYPosition(:,2),1,colourVector(:,1))
        hold on
    else
        title({['2-D Plot of Particle Trajectories'],['Seth Thompson | 101031310']})
        scatter(particleXPosition(:,2),particleYPosition(:,2),1,colourVector(:,1))
        hold off
    end
    
    pause(0.001)
    
    % Updating the first column of each position vector so the position of
    % each particle can be incremented
    particleXPosition(:,1) = particleXPosition(:,2);
    particleYPosition(:,1) = particleYPosition(:,2);
    
    counter = counter + 1;
end

% 1D) What is the relationship between the electron drift current density
% and average carrier velocity? Give the formula for current and generate a
% plot of the current overtimeinthe x direction. Assume an electron
% concentration of 10^15 1/cm^2. Note that the electron concentration is used
% to find out how many electrons each of the model’s “particles” represent.
% Comment on the current behaviour over time. 

% It is known that the equation for drift current density is...
% J = N*(charge/(density))*mean(Vn)
% From there, the drift current density can be calculated from the average velocity at
% a given point in time. The code made to do this is shown next.
concentration = 10^19; %1/m^2

% Simulate a larger amount of electrons here for more accurate results.
elecAmount2 = 10000;

% Initialize random velocities and displacements once again for the new simulation.
randXV = randn(elecAmount2,2)*sqrt((kb*Temperature)/mEffective);
randYV = randn(elecAmount2,2)*sqrt((kb*Temperature)/mEffective);

% Assigning a particle a random position on the XY plane within the maximum
% and minimum width and length.
particleXPosition = regionLength.*rand(elecAmount2,2);
particleYPosition = regionWidth.*rand(elecAmount2,2);

% First column is the initial position, second column is the next position,
% think xn and xn+1
particleXPosition(:,1) = particleXPosition(:,2);
particleYPosition(:,1) = particleYPosition(:,2);

% Create a loop that will calculate a new velocity for each particle after
% each time-step and also calculate the average velocity at that point.
% Creating a loop that will calculate each particles position as time goes
% on. Also making a counter that will keep track of the amount of
% iterations done. 
counter = 0;

% Calculating the displacement of the electron for each time step
particleXDisplacement = timeStep*randXV(:,1);
particleYDisplacement = timeStep*randYV(:,1);

% creating variables that will hold the average velocities in the x and y
% direction
Vx = 0;

for m = 1:1000

    randXV = randXV + acceleration*timeStep;
    
    % Calculating the average velocity of all particles in the x direction
    Vx(m) = mean(randXV(:,2));
    % Updating the position of the particle
    for n = 1:elecAmount2
        % First checking to see if the electron scatters, if it does then a
        % new velocity must be generated and from there a new displacement
        if (rand < PScat)
            thetaScatter = (2*pi).*rand(1);
            randXV(n,:) = randn(1)*sqrt((kb*Temperature)/mEffective) + acceleration*timeStep;
            randYV(n,:) = randn(1)*sqrt((kb*Temperature)/mEffective);
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleYDisplacement(n,:) = timeStep*randYV(n,1);           
            
            % Saving the new positions before the next collision
            beforeCollideX = particleXPosition(:,1);
            beforeCollideY = particleYPosition(:,1);
        end
        % Checking the boundary conditions for the left and right sides of
        % the plot. If the particle is going to pass a border on the left
        % or right side then simply move it to the left or right.
        if (particleXPosition(n,1) + particleXDisplacement(n) > 2e-7)
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1) - 2e-7;
        elseif (particleXPosition(n,1) + particleXDisplacement(n) < 0)
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1) + 2e-7;
        else
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1);
        end
        
        % Checking boundary conditions for the top and bottom of the plot.
        % If the particle is about to cross the top of the plot, flip the
        % y, displacement value and have it move in the opposite direction;
        % the opposite applies for the bottom of the plot.
        if ((particleYPosition(n,1) + particleYDisplacement(n) > 1e-7) || (particleYPosition(n,1) + particleYDisplacement(n) < 0))
            particleYDisplacement(n) = -particleYDisplacement(n);
            particleYPosition(n,2) = particleYDisplacement(n) + particleYPosition(n,1);
        else
            particleYPosition(n,2) = particleYDisplacement(n) + particleYPosition(n,1);
        end
    end
    timeVector(m) = timeStep*m;
end

% Creating a plot of the current density over time.
currentDen = (regionWidth*(eCharge)*(concentration)).*Vx;
figure(4)
plot(timeVector,currentDen)
title({'Current Density in the X-Direction','Seth Thompson | 101031310'})
xlabel('Time (s)')
ylabel('Current Densiy (A/m)')
grid on

% Looking at the plot of the current density, the results indicate that it
% is converging to a certain value as time goes on, which makes sense,
% The value it is converging to seems reasonable as well.

% 1E) As before generate the density and temperature maps at the end of the
% simulation.

% To do this, the hist3 functions will be used to create a 3-d histogram
% with the bivariate data (X and Y positions of each particle).

% Making the electron density plot.

figure(5)
hist3([particleXPosition(:,1),particleYPosition(:,1)],[10,20])
title({['Electron Density Map'],['Seth Thompson | 101031310']})
xlabel('X-Axis Segments')
ylabel('Y-Axis Segments')
zlabel('Number of Particles')

% Calculating the final velocities of eqch electron
randV = (randXV(:,1).^2 + randYV(:,1).^2).^(0.5);

% Code to make a temperature map with colours.
% To start, the kinetic energy of each particle at the end of the
% simulation must be calculated.

kineticEtempPlot = (mEffective.*(randV.^2))/2;
individualTemp = kineticEtempPlot./kb;
individualTempC = individualTemp - 273;

% Scaling xposition and y position vectors up so they can be used to access
% difierent parts of the array defined previously.

XPosScaled = particleXPosition(:,1).*1e9;
YPosScaled = particleYPosition(:,1).*1e9;

% Rounding each scaled vector so its elements can be used to access parts
% of an array
XPosScaled = round(XPosScaled,-1)./10;
YPosScaled = round(YPosScaled,-1)./10;

% A 10x20 matrix must now be made holding the average of each temperature
% in a 'bin'.

temperatureArray = zeros(round(regionWidth/1e-9)/10,round(regionLength/1e-9)/10);

% Creating a loop that will store the previously defined array with the sum
% of the temperature values in the array.

for r = 1:electronAmount
    XPosScaled = round(particleYPosition(r,1).*1e9,-1)/10;
    YPosScaled = round(particleXPosition(r,1).*1e9,-1)/10;
    
    if(XPosScaled <= 0 || YPosScaled <= 0)
        XPosScaled = 1;
        YPosScaled = 1;
    end
    temperatureArray(XPosScaled,YPosScaled) = temperatureArray(XPosScaled,YPosScaled) + individualTempC(r);
end

% Storing the number of elements in each bun in a matrix.
divisionElements = hist3([particleXPosition(:,1),particleYPosition(:,1)],[10,20]);

% Creating another loop that will divide the temperatures in each element
% by the number of particles there.

for p = 1:round(regionWidth/1e-9)/10
    for q = 1:round(regionLength/1e-9)/10
        if(divisionElements(p,q) == 0)
            temperatureArray(p,q) = 0;
        else
            temperatureArray(p,q) = temperatureArray(p,q)/divisionElements(p,q);
        end
    end
end

% Plotting the histogram for the temperature
figure(6)
bar3(temperatureArray)
title({['Temperature Map'],['Seth Thompson | 101031310']})
xlabel('X-Axis Segments')
ylabel('Y-Axis Segments')
zlabel('Temperature (C)')