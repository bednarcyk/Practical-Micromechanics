% -- Written by Mathew Schey, U Mass Lowell student, Dec. 2021
sz = get(0, 'screensize');
h = figure("position",[sz(1)+100, sz(2)+200, 0.9*sz(3)+100, sz(4)/2+200]);
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'failure.gif';
a1 = dir('Damage Snapshot*.bmp');
[xx, order] = sort([a1(:).datenum],'ascend');
a = a1(order);
for n = 1:size(a,1)
    
    I = imread(a(n).name);
    imshow(I)
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256); % -- Matlab
    % [imind,cm] = rgb2ind(im); % -- Octave
    
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end
    
end
close(h);
fclose('all');
