% Read the input image
image = imread('test.jpg'); % Replace 'your_image.jpg' with your image file name

% Convert the image to grayscale
grayImage = rgb2gray(image);
ave = sum(grayImage(:))/numel(grayImage);
grayImage = grayImage-ave;
grayImage(grayImage < 20) = 0;

% Find connected components (regions) in the matrix
CC = bwconncomp(grayImage);

area = zeros(1,CC.NumObjects);
origamis = zeros(1,5);
maximum = max(grayImage(:));


% Initialize a figure
figure;
alphaValue = 0.4;
imshow(grayImage,[]);
hold on;
% Loop through the connected components
for i = 1:CC.NumObjects
    % Get the indices (linear indices) of pixels belonging to the current region
    regionIndices = CC.PixelIdxList{i};
    area(1,i) = length(regionIndices);  

    % Extract intensity values for the current region
    regionIntensityValues = grayImage(regionIndices);
    maximum_inregion = max(regionIntensityValues);

    [row, col] = ind2sub(size(grayImage), regionIndices(1));

    first = 100;
    second = 210;
    thrid = 280;
    fourth = 400;
    
    DLintensity = maximum*0.7;
    if area(1,i)>first && area(1,i)<=second
        scatter(col, row, 50, 'filled', 'MarkerFaceColor', 'r','MarkerFaceAlpha', alphaValue);
        origamis(1,1)=origamis(1,1)+1;
    elseif area(1,i)>second && area(1,i) <= thrid
        if maximum_inregion <= DLintensity
            scatter(col, row, 50, 'filled', 'MarkerFaceColor', 'g','MarkerFaceAlpha', alphaValue);
            origamis(1,2)=origamis(1,2)+1;
        elseif maximum_inregion > DLintensity
            scatter(col, row,50, 'filled', 'MarkerFaceColor', 'y','MarkerFaceAlpha', alphaValue);
            origamis(1,3)=origamis(1,3)+1;
        end
    elseif area(1,i)>thrid && area(1,i) <=460 && maximum_inregion <= DLintensity
        scatter(col, row, 50, 'filled', 'MarkerFaceColor', 'b','MarkerFaceAlpha', alphaValue);
        origamis(1,4)=origamis(1,4)+1;   
    end
end

hold off;
disp(origamis);

%histogram(area, 'BinWidth', 3, 'FaceColor', 'b', 'EdgeColor', 'k');