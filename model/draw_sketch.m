% Takes in output from pencil_sketch and combines all ridges into image

X = [];
Y = [];
files = ls;

for i = 1:length(ls)
    if contains(files(i,:), 'A_')
        load(files(i,:))
        X = [X Xs];
        Y = [Y Ys];
    end
end

% Produce sketch
sketch = zeros(360,360);
sketch(sub2ind(size(sketch),Y,X)) = 1;

figure()
imagesc(sketch)
set(gca,'YDir','normal')
axis image

% Widen sketch
sketch2 = imgaussfilt(sketch,1);
threshold = 0.2;
sketch2(sketch2 > threshold) = 1;
sketch2(sketch2 <= threshold) = 0;
sketch2 = imgaussfilt(sketch2,5);

figure()
imagesc(sketch2)
set(gca,'YDir','normal')
axis image
