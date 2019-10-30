% Creates a pencil sketch by interpolating between points you pick.
% Automatically discards first click when zooming. 
% Need draw_sketch to create an image. 

path = 'Your path';
image = fitsread([path  'IMAGE_141_J8.9.fits']);
image = image/max(max(image));

high_pass_data = image - imgaussfilt(image,5);
high_pass_data = min(max(high_pass_data,0),0.012);
high_pass_data = imgaussfilt(high_pass_data, 1);

figure()
imagesc(high_pass_data)
set(gca,'YDir','normal')
axis image

Xs = [];
Ys = [];
done = "";

while done ~= "y"
    [x, y] = getpts;
    accept = input('Accept? (y/n) ','s');
    if accept == 'y'
        knots = 1:length(x(2:end))+3;
        curve = spmak( knots, [x(2:end)'; y(2:end)'] );
        t = linspace(3,length(x(2:end))+1,length(x(2:end))*100);
        values = fnval(curve,t);
        
        Xs = [Xs, values(1,:)];
        Ys = [Ys, values(2,:)];
    end
    done = input('Done? (y/n) ','s');
end

% Produce sketch
Xs = round(Xs);
Ys = round(Ys);

sketch = zeros(size(image));
sketch(sub2ind(size(sketch),Ys,Xs)) = 1;

figure()
imagesc(sketch)
set(gca,'YDir','normal')
axis image

% Widen sketch
sketch2 = imgaussfilt(sketch,1);
threshold = 0.2;
sketch2(sketch2 > threshold) = 1;
sketch2(sketch2 <= threshold) = 0;

figure()
imagesc(sketch2)
set(gca,'YDir','normal')
axis image

name = input('Name spline: ','s');
if name ~= "no"
    save("Your path"+name+".mat",'Xs','Ys')
end

close all