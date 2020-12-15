degree = 30;
init = [1 0 0]';
step = 360/degree;

x = zeros(1, step);
y = zeros(1, step);
z = zeros(1, step);
x(1) = init(1);
y(1) = init(2);
z(1) = init(3);

rotationMatrix = [cosd(degree), -sind(degree), 0; sind(degree), cosd(degree), 0; 0 0 1];
scatter3(x(1), y(1), z(1));
hold on;
for numIdx = 2: step
    point = [x(numIdx - 1), y(numIdx - 1), z(numIdx - 1)]';
    updatedPoint = rotationMatrix*point;
    x(numIdx) = updatedPoint(1);
    y(numIdx) = updatedPoint(2);
    z(numIdx) = updatedPoint(3);
    scatter3(x(numIdx), y(numIdx), z(numIdx));
end
hold off;