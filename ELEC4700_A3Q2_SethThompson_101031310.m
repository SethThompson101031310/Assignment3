%% ELEC 4700 Assignment 3 Question Two Work

% Name: Seth Thompson
% Student Number: 101031310

close all
clear
clc

% Starting this code by putting in my solution from assignment 2.

% Defining the length and width of the box along with Vo

Length = 2;
Width = 1;
Vo = 1;
% Defining the values for sigma both inside and outside of the boxes.

sigmaIN = 1e-2;
sigmaOUT = 1;

% Defining the dimensions of each of the boxes (wb and lb in figure 3)

WidthBox = 0.4;
LengthBox = 0.4;

% Defining the number of elements that will be in each part of the matrix

nx = 100*Length;
ny = 100*Width;

% Defining the conductivity map

conductivity = zeros(ny,nx);

for k = 1:ny
    for l = 1:nx
        % If the element being accesed in the conductivity matrix
        % IS within one of the boxes, set its value to the lower
        % conductivity value
        if(l >= nx*WidthBox && l <= nx-nx*WidthBox && (k >= ny-ny*LengthBox || k <= ny*LengthBox))
            conductivity(k,l) = sigmaIN;
        % Else, put in the higher value
        else
            conductivity(k,l) = sigmaOUT;
        end
    end
end

G = sparse(nx*ny,nx*ny);
B = zeros(nx*ny,1);

% Populating the G matrix
for l = 1:nx
    for k = 1:ny
        
        % Node mapping to put entries into the correct place
        n = k + (l - 1)*ny;
        
        % Calculating deltas in the x and y direction
        nxm = k + (l - 2)*ny;
        nxp = k + l*ny;
        nym = (k - 1) + (l - 1)*ny;
        nyp = (k + 1) + (l - 1)*ny;
        
        % Begin inputting all of the correct entires!
        if(l == 1 || l == nx) % Left or right side of the region, set entries to 1
            G(n,n) = 1;
        elseif (k == 1) % We are along the bottom of the region apply finite difference method as needed
            
            entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
            entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
            entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;
            
            G(n,n) = -(entryYup + entryXup + entryXdown);
            G(n,nyp) = entryYup;
            G(n,nxp) = entryXup;
            G(n,nxm) = entryXdown;
            
        elseif (k == ny) % We are along the top of the region
            
            entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
            entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
            entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;
            
            G(n,n) = -(entryYdown + entryXup + entryXdown);
            G(n,nym) = entryYdown;
            G(n,nxp) = entryXup;
            G(n,nxm) = entryXdown;
        else % else, apply finite differnce as needed without worrying about going out of bounds...
            
            % Storing elements from conductivity matrix that will be mapped
            % into the G matrix
            entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
            entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
            entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
            entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;
            
            G(n,n) = -(entryYup + entryYdown + entryXup + entryXdown);
            G(n,nyp) = entryYup;
            G(n,nym) = entryYdown;
            G(n,nxp) = entryXup;
            G(n,nxm) = entryXdown;
            
        end
    end
end

% Populating the B vector next...
for l = 1:nx
    for k = 1:ny
        % Node mapping to put entries into the correct place
        n = k + (l - 1)*ny;
        
        % Are we along the left side? if so set the value to Vo
        if (l == 1) 
            B(n) = Vo;
        end
        
        % Anywhere also it should be zero, but that was defined by using
        % the zeros function to make the vector.
    end
end

% Obtaining the solution, V, from V = G\B

V = G\B;

% Moing the Solution V from the signle vector to a matrix 

for l = 1:nx
    for k = 1:ny
        % Node mapping to put entries into the correct place
        n = k + (l - 1)*ny;
        
        MatrixV(k,l) = V(n);
    end
end

% plotting the voltage with surf
figure(1)
surf(MatrixV)
xlabel('Length')
ylabel('Width')
zlabel('Voltage (V)')
title({'Voltage Over the Region','Seth Thompson | 101031310'})
% Pausing so the plot will plot before the rest of the code runs
pause (1)
% The simulated voltage plot makes sense, it slowly drops before getting
% near the contacts and the drops rapidly once it passes through them, then
% goes back to dropping slowly to 0.

% Using the Gradient Function to obtain the electric field from the voltage
% (E = -del(V))

[Ex,Ey] = gradient(-MatrixV,1e-9);

% Using the Quiver Function to plot the electric field over the region
figure(2)
quiver(Ex,Ey)
xlabel('Length')
ylabel('Width')
title({'Electric Field over The Region','Seth Thompson | 101031310'})
ylim([0,100])
xlim([0,200])
% Pausing so the plot will plot before the rest of the code runs
pause (1)
% The simulated electric field plot makes sense, it is in the direction of
% higher to lower potential.

% The electric field matrix calculated before will be used to apply the
% correct amount of acceleration to each electron being simulated. The code
% from part one of this assignment will be re-used here.

eCharge = 1.602e-19; % Columbs
mRest = 9.109e-31; % kilograms
mEffective = 0.26*mRest; % kilograms

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

electronAmount = 10000; % Amount of electrons being simulated

% Assigning a particle a random position on the XY plane within the maximum
% and minimum width and length.
particleXPosition = regionLength.*rand(electronAmount,2);
particleYPosition = regionWidth.*rand(electronAmount,2);

for y = 1:electronAmount
    
    % Has a particle been defined to be inside one of the boxes? then
    % redefine its starting position randomly until its not
    while ((particleXPosition(y,2) > 0.8e-7 && particleXPosition(y,2) < 1.2e-7) && ...
            (particleYPosition(y,2) > 0.6e-7 ||  particleYPosition(y,2) < 0.4e-7))
        particleXPosition(y,2) = regionLength.*rand(1,1);
        particleYPosition(y,2) = regionWidth.*rand(1,1);
    end

end

% First column is the initial position, second column is the next position,
% think xn and xn+1
particleXPosition(:,1) = particleXPosition(:,2);
particleYPosition(:,1) = particleYPosition(:,2);

% Timestep calculated from sqrt(100nm^2 + 200nm^2)/1000, which is approx
% 0.22nm, and then 0.22nm/thermalVel
timeStep = 0.22e-9/thermalVel; % seconds

% for this assignemtn, make time step a little smaller
timeStep = timeStep/2;

% Assigning a particle a random velocity using the randn function
% (normally distributed random values).
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

for m = 1:1000
    
    % Updating the position of the particle
    for n = 1:electronAmount
        
        % Must update velocities based on the electric field stregth where
        % the particle is. This is done next.
        if(round(particleYPosition(n,1)*1e9) <= 0 || round(1e9*particleXPosition(n,1)) <= 0)
            Forcex = Ex(ceil(particleYPosition(n,1)*1e9),ceil(1e9*particleXPosition(n,1))).*eCharge; % Newtons
            Forcey = Ey(ceil(particleYPosition(n,1)*1e9),ceil(1e9*particleXPosition(n,1))).*eCharge; % Newtons
        elseif(round(particleYPosition(n,1)*1e9) > 100 || round(1e9*particleXPosition(n,1)) > 200)
            Forcex = Ex(floor(particleYPosition(n,1)*1e9),floor(1e9*particleXPosition(n,1))).*eCharge; % Newtons
            Forcey = Ey(floor(particleYPosition(n,1)*1e9),floor(1e9*particleXPosition(n,1))).*eCharge; % Newtons
        else
            Forcex = Ex(round(particleYPosition(n,1)*1e9),round(1e9*particleXPosition(n,1))).*eCharge; % Newtons
            Forcey = Ey(round(particleYPosition(n,1)*1e9),round(1e9*particleXPosition(n,1))).*eCharge; % Newtons
        end 
        accelerationx = Forcex/mEffective;
        accelerationy = Forcey/mEffective;
        randXV(n,1) = randXV(n,1) + accelerationx*timeStep;
        randYV(n,1) = randYV(n,1) + accelerationy*timeStep;
        
        % Checking to see if the electron scatters, if it does then a
        % new velocity must be generated and from there a new displacement
        if (rand < PScat)
            randXV(n,:) = randn(1)*sqrt((kb*Temperature)/mEffective);
            randYV(n,:) = randn(1)*sqrt((kb*Temperature)/mEffective);
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleYDisplacement(n,:) = timeStep*randYV(n,1);           
            
            % Saving the new positions before the next collision
            beforeCollideX = particleXPosition(:,1);
            beforeCollideY = particleYPosition(:,1);
        else
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleYDisplacement(n,:) = timeStep*randYV(n,1);
        end
        % Checking the boundary conditions for the left and right sides of
        % the plot. If the particle is going to pass a border on the left
        % or right side then simply move it to the left or right.
        if (particleXPosition(n,1) + particleXDisplacement(n) > 2e-7)
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1) - 2e-7;
        elseif (particleXPosition(n,1) + particleXDisplacement(n) < 0)
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1) + 2e-7;
        else
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1);
        end
        
        % Checking boundary conditions for the top and bottom of the plot.
        % If the particle is about to cross the top of the plot, flip the
        % y, displacement value and have it move in the opposite direction;
        % the opposite applies for the bottom of the plot.
        
        % For part three of the assignment, more boundary conditions will
        % be added here so that the particles bounce off of the boxes in
        % the simulation.
        % Does it hit the top or bottom?
        if ((particleYPosition(n,1) + particleYDisplacement(n) > 1e-7) || (particleYPosition(n,1) + particleYDisplacement(n) < 0))
            particleYDisplacement(n) = -particleYDisplacement(n);
            particleYPosition(n,2) = particleYDisplacement(n) + particleYPosition(n,1);
        % ... right side of the bottom box?
        elseif (particleXPosition(n,1) > 0 && particleXPosition(n,1) < 0.8e-7 && ...
                (particleXPosition(n,1) + particleXDisplacement(n) > 0.8e-7) && ...
                particleYPosition(n,1)> 0 && particleYPosition(n,1) < 0.4e-7)
            particleXDisplacement(n) = -particleXDisplacement(n);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1);
        % ... right side of the top box?
        elseif (particleXPosition(n,1) > 0 && particleXPosition(n,1) < 0.8e-7 && ...
                particleYPosition(n,1) > 0.6e-7 && particleYPosition(n,1) < 1e-7 && ...
                (particleXPosition(n,1) + particleXDisplacement(n) > 0.8e-7))
            particleXDisplacement(n) = -particleXDisplacement(n);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1);
        % ... left side of the top box?
        elseif(particleXPosition(n,1) > 1.2e-7 && particleXPosition(n,1) < 2e-7 && ...
                particleYPosition(n,1) > 0.6e-7 && particleYPosition(n,1) < 1e-7 && ...
               (particleXPosition(n,1) + particleXDisplacement(n) < 1.2e-7))
            particleXDisplacement(n) = -particleXDisplacement(n);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1);
        % ... left side of the bottom box?
        elseif(particleXPosition(n,1) > 1.2e-7 && particleXPosition(n,1) < 2e-7 && ...
                particleYPosition(n,1) > 0 && particleYPosition(n,1) < 0.4e-7 && ...
               (particleXPosition(n,1) + particleXDisplacement(n) < 1.2e-7))
            particleXDisplacement(n) = -particleXDisplacement(n);
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1);
        % ... bottom of the top box?
        elseif(particleXPosition(n,1) > 0.8e-7 && particleXPosition(n,1) < 1.2e-7 && ...
                particleYPosition(n,1) > 0.4e-7 && particleYPosition(n,1) < 0.6e-7 && ...
               (particleYPosition(n,1) + particleYDisplacement(n) < 0.4e-7))
            particleYDisplacement(n) = -particleYDisplacement(n);
            particleYPosition(n,2) = particleYDisplacement(n) + particleYPosition(n,1); 
        % ... top of the bottom box?
        elseif(particleXPosition(n,1) > 0.8e-7 && particleXPosition(n,1) < 1.2e-7 && ...
                particleYPosition(n,1) > 0.4e-7 && particleYPosition(n,1) < 0.6e-7 && ...
               (particleYPosition(n,1) + particleYDisplacement(n) > 0.6e-7))
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
        scatter(particleXPosition([1:10],2),particleYPosition([1:10],2),1,colourVector([1:10],1))
        hold on
        % Plotting the bottom box
        plot([0.8e-7,0.8e-7],[0,0.4e-7],'k',[0.8e-7,1.2e-7],[0.4e-7,0.4e-7],'k',[1.2e-7,1.2e-7],[0.4e-7,0],'k',[0.8e-7,1.2e-7],[0,0],'k')
        
        % Plotting the top box.
        plot([0.8e-7,0.8e-7],[1e-7,0.6e-7],'k',[0.8e-7,1.2e-7],[0.6e-7,0.6e-7],'k',[1.2e-7,1.2e-7],[0.6e-7,1e-7],'k',[0.8e-7,1.2e-7],[1e-7,1e-7],'k')
        
        % Labelling axis and defining limtis as needed
        title({['2-D Plot of Particle Trajectories'],['Seth Thompson | 101031310']})
        xlabel('X-Axis (m)')
        ylabel('Y-Axis (m)')
        xlim([0 200e-9])
        ylim([0 100e-9])
    elseif (counter > 0 && counter < 1000)
        % The title of the plot will update after each iteration with the
        % new average for the semiconductors temperature
        title({['2-D Plot of Particle Trajectories'],['Seth Thompson | 101031310']})
        scatter(particleXPosition([1:10],2),particleYPosition([1:10],2),1,colourVector([1:10],1))
        hold on
    else
        title({['2-D Plot of Particle Trajectories'],['Seth Thompson | 101031310']})
        scatter(particleXPosition([1:10],2),particleYPosition([1:10],2),1,colourVector([1:10],1))
        hold off
    end
    
    pause(0.001)
    
    % Updating the first column of each position vector so the position of
    % each particle can be incremented
    particleXPosition(:,1) = particleXPosition(:,2);
    particleYPosition(:,1) = particleYPosition(:,2);
    
    counter = counter + 1;
end

% As can be seen from the simulation results, it works as intended except
% for the odd case the a particle comes towards a box at its corners, and
% it clips through.
