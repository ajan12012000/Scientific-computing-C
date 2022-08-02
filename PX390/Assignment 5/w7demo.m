warning off; close all; clear all;
% % Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman');
set(0,'DefaultAxesFontSize', 20);
% % Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20);
%% Numerical solution
% Load the data and the simulations parameters
data0=load('heq2d.dat');
sim_params=load('heq2d_params');
nx = sim_params(1);
% The data can be reshaped to a more convinient format
data=reshape(data0, size(data0,1)/nx, nx, size(data0,2));
% We can visualise the time evolution of temperature by plotting the
% temperature as a color along the grid and at increasing time
figure;
pcolor(data(:,:,3));
shading flat; colormap(jet(100)); newcolmap=colormap;
colorbar('EastOutside');
xlabel('Position x'); ylabel('Position y');
%% And now a linear plot of a profile in a cut through the domain in one
%  particular direction, here cut at x=nx/2 + 1
x0 = nx/2 + 1;
figure;
plot(data(x0,:,3));
xlabel("Position y"); ylabel("Temperature");