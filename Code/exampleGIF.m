% Create the figure and axes
figure;
histogram_axes = axes;

% Set up the animated histogram
animated_hist1 = animatedline(histogram_axes, 'Color', 'r');
animated_hist2 = animatedline(histogram_axes, 'Color', 'g');
animated_hist3 = animatedline(histogram_axes, 'Color', 'b');
animated_hist4 = animatedline(histogram_axes, 'Color', 'r');
animated_hist5 = animatedline(histogram_axes, 'Color', 'g');
animated_hist6 = animatedline(histogram_axes, 'Color', 'b');


xlabel('types');
ylabel('intensity');
title('origami types distribution across time');

% Loop through the time steps and update the histogram
for t = 1:length(time)
    % Update the data for the histograms
    histogram(listm(1:t), 'FaceColor', 'r', 'EdgeColor', 'none', 'Parent', histogram_axes);
    histogram(listd(1:t), 'FaceColor', 'g', 'EdgeColor', 'none', 'Parent', histogram_axes);
    histogram(listtrimers(1:t), 'FaceColor', 'g', 'EdgeColor', 'none', 'Parent', histogram_axes);
    histogram(listte(1:t), 'FaceColor', 'g', 'EdgeColor', 'none', 'Parent', histogram_axes);
    histogram(listp(1:t), 'FaceColor', 'g', 'EdgeColor', 'none', 'Parent', histogram_axes);
    histogram(listhexmers(1:t), 'FaceColor', 'g', 'EdgeColor', 'none', 'Parent', histogram_axes);
    
    % Capture the frame
    frame = getframe(gcf);
    image_data = frame2im(frame);
    
    % Convert the frame to indexed GIF format
    [imind, cm] = rgb2ind(image_data, 256);
    
    % Write the frame to the GIF file
    if t == 1
        imwrite(imind, cm, 'histogram.gif', 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
    else
        imwrite(imind, cm, 'histogram.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end
    
    % Clear the axes for the next frame
    cla(histogram_axes);
end