function [A, S, C] = preprocess_image(image_path, nbox, f)
    % PREPROCESS_IMAGE  Convert an image into initial density grids for
    % the tip/stalk/collision model, and plot the three densities.
    %
    % Inputs:
    %   image_path - path to the source image (e.g. 's3.png')
    %   nbox       - number of boxes to divide the lateral field into
    %                (spatial grid resolution)
    %   f          - accepted for interface compatibility with the
    %                caller, but not used inside this function
    %
    % Outputs:
    %   A, S, C    - active-tip, stalk, and inactive/collided density
    %                grids (nbox x nbox)

    % Read the image
    image = imread(image_path);

    % Classify each pixel by color into the tracked populations:
    %   red    -> active tips
    %   black  -> stalks
    %   yellow -> inactive/collided entities
    % (green is also computed but not used further)
    black_channel = sum(image, 3) == 0; % black assumed to be RGB [0 0 0]
    red_channel = image(:,:,1) - image(:,:,2) == 255;
    yellow_channel = image(:,:,2) - image(:,:,3) == 255;
    green_channel = image(:,:,2) - image(:,:,1) > 50;

    % Initial variable setup
    I = size(image);
    npix = fix(I(1) / nbox); % pixels per box; assumes the image is square, fix() keeps the integer part

    % Initialize matrices for storing density values
    A = zeros(nbox,nbox);
    S = zeros(nbox,nbox);
    C = zeros(nbox,nbox);

    % Calculate the density in each box by counting matching pixels
    % (indexing A(j,i,1) below is equivalent to A(j,i) since A is 2-D;
    % MATLAB allows trailing singleton indices)
    for j = 1:nbox
        nr = npix * (j - 1);
        for i = 1:nbox
            nc = npix * (i - 1);
            for k = 1:npix
                for m = 1:npix
                    A(j, i) = A(j, i, 1) + red_channel(nr + k, nc + m);
                    S(j, i) = S(j, i, 1) + black_channel(nr + k, nc + m);
                    C(j, i) = C(j, i, 1) + yellow_channel(nr + k, nc + m);
                end
            end
        end
    end
    X = linspace(0,1,nbox);% row vector of nbox evenly spaced points between 0 and 1
    Y = linspace(0,1,nbox);% in this model's scale, 0->1 corresponds to about 0->1 mm
    dx = X(2) - X(1);

    % Normalize pixel counts into densities: each tip corresponds to 953 pixels
    A = A/(953*dx*dx);
    S = S/(953*dx*dx);
    C = C/(953*dx*dx);

    %% Plot the three initial densities (figures 1-3)
    figure(1)
    surf(X,Y,A(:,:)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar

    figure(2)
    surf(X,Y,S(:,:)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar

    figure(3)
    surf(X,Y,C(:,:)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar
end
