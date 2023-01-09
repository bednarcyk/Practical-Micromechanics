% -- Written by Mathew Schey, U Mass Lowell student, Dec. 2021
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'failure.gif';
a = dir('Damage Snapshot*.bmp');
a = struct2table(a);
a = sortrows(a,'date');
a = table2struct(a);
for n = 1:size(a,1)
    I = imread(a(n).name);
    imshow(I) 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); s
    % Write to the GIF File 
    if n == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1); 
    end 
end
 
fclose('all');
