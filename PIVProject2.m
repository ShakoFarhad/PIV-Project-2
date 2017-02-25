addpath(genpath('/mn/kadingir/domt_143903/HydrolabPIV/src'));
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/NewCoordinates/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run1(0.05)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run2(0.05)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run3(0.05)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run1(0.1)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run2(0.1)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run3(0.1)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run1(0.2)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run2(0.2)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run3(0.2)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run1(0.3)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run2(0.3)/PIV-Images');
addpath('/mn/marduk/data/MEK4600_2017/gruppe1/Data21022017/Run3(0.3)/PIV-Images');
javaaddpath('/mn/kadingir/domt_143903/HydrolabPIV/src/measures');
javaaddpath('/mn/kadingir/domt_143903/HydrolabPIV/src/interp');


% Read in wave images
im1 = imread('Crest_03_1A.tif');
mask1 = imread('MaskCrest_03_1A.png');
im2 = imread('Crest_03_1B.tif');
mask2 = imread('MaskCrest_03_1B.png');
coord = imread('coord.tif');

imagesc(coord) % Takes the picture with the marked coordinate points and shows it.
% Select reference points in pixel coordinate
h=impoly; % Creates lines between points.
pixel = h.getPosition; % Gets the position of the lines' points.

% Refine pixel positions (optional)
bw = im2bw(coord); % makes the coordinate image black and white
cc = bwconncomp(bw); % 
stats = regionprops(cc,'Centroid');
xc = vertcat(stats.Centroid);
idx = knnsearch(xc,pixel);
pixel = xc(idx,:);

% Define matching reference points in world coordinate
[wx, wy] = ndgrid((0:1:10)*0.02, (-1:-1:-9)*0.02);
world = [wx(:) wy(:)];

% Create coordinate transformation
[tform,err,errinv] = createcoordsystem(pixel,world,'linear');

g = 9.81; % m/s^2
f = 1.425; % Frequency 1/s
a = 0.0147; % m
w = 2*pi*f; % 1/s
k = 8.172; % 1/m
z = linspace(-0.2, 0,200);
U_ex = a*w*exp(k*z); %Stokes
%Linear


for n = [72, 64, 60, 56, 52, 44] % Subwindow size in pixels
    figure
    counter = 0;
    for k = [2, 3, 4, 5]
        overlap = 0.5; % Percentage overlap 
        m = floor(n/k); %Search area size in pixels
        
        % First pass with masks
        opt1 = setpivopt('range',[-m m -m m],'subwindow',n,n,overlap);
        piv1 = normalpass([],im1,mask1,im2,mask2,opt1);
        dt = 0.01;
        [U1,V1,x,y] = replaceoutliers(piv1,mask1&mask2);
        [Uw,Vw,xw,yw] = pixel2world(tform,U1,V1,x,y,dt);
   
        % Finding the lowest vertical velocity coordinate, ie. finding the
        % index for the column of piv1.V, where piv1.V has the lowest velocity.
        columnSize = 0;
        minimumColumnSize = 100000;
        columnIndex = 0;
        for i = 1:size(piv1.V,2)
            columnSize = sum(abs(piv1.V(:,i)))/size(piv1.V,2);
            if(columnSize < minimumColumnSize)
                minimumColumnSize = columnSize;
                columnIndex = i;
            end
        end
    
        % Plotting the different U velocities in world frame in each subplot.
        counter = counter + 1;
        subplot(2,2,counter)
        plot(Uw(:,columnIndex),yw(:,columnIndex))
        hold on
    
        % Plotting the analytical solution
        plot(U_ex, z)
        hold on

        title(['k*a = ', num2str(k*a), '. Subwindow size = ', num2str(n), ', search range = ', num2str(m), ' and overlap = ', num2str(overlap*100), '%.'])
        legend('PIV U', 'U exact')
        
        xlabel(' u [ m/s ] ');
        ylabel(' z [ m ] ');
    end
end

[U2,V2,x2,y2] = replaceoutliers(piv1);

figure
load wind
cav = curl(x2,y2,U2,V2); %plotting vorticity
pcolor(x2,y2,cav); shading interp
hold on;
quiver(x2,y2,U2,V2, 3, 'y')
hold off
colormap winter
xlabel(' u ');
c = colorbar;
ylabel(c,' z ');
axis ij

[U2,V2,x2,y2] = replaceoutliers(piv1,mask1&mask2);
figure
quiver(x2,y2,U2,V2, 3)
xlabel(' u ');
ylabel(' z ');
axis ij