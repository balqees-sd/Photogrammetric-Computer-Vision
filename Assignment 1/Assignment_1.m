% Task 3

% ==============================================================================

% Part a

% Points represented in homogenous coordinate system
x = [2; 3; 1];
y = [-4; 5; 1];

% The cross product of x and y will give us the connecting line between them
l = cross(x, y);
disp('Connecting line between the two points:'), disp(l);

% ==============================================================================

% Part b

% The translation matrix mT for movement in the direction (6, -7)
mT = [
1 0 6;
0 1 -7;
0 0 1
];

% Multiplying the two points by the translation matrix:
x1 = mT * x;
y1 = mT * y;

disp('The points x and y after translation are:');
disp('x1:'), disp(x1);
disp('y1:'), disp(y1);

% The rotation matrix mR for rotation by 15 degrees
phi = deg2rad(15);
mR = [
cos(phi) -sin(phi) 0;
sin(phi) cos(phi) 0;
0 0 1;
];

% Multiplying the points x1 and y1 by the rotation matrix:
x2 = mR * x1;
y2 = mR * y1;

disp('The points x1 and y2 after rotation are:');
disp('x2:'), disp(x2);
disp('y2:'), disp(y2);

% The scale matrix mS for scaling by 8 units
lambda = 8;
mS = [
8 0 0;
0 8 0;
0 0 1;
];

% Multiplying the points x1 and y1 by the rotation matrix:
x3 = mS * x2;
y3 = mS * y2;

disp('The points x2 and y3 after scaling are:');
disp('x3:'), disp(x3);
disp('y3:'), disp(y3);

% ==============================================================================

% Part c

% Performing transformations with one matrix
mTransformation = mS * mR * mT;

% Taking the inverse and then the transpose of the transformation matrix
mTransformationInv = inv(mTransformation);
mTransformationInvTrans = transpose(mTransformationInv);

% Multiplying mTransformationInvTrans with the line
lTransformation = mTransformationInvTrans * l;

xIntercept = (-lTransformation(3)/lTransformation(1));
yIntercept = (-lTransformation(3)/lTransformation(2));

disp('The transformation applied to the line results in:');
disp(lTransformation);

disp('X intercept of the line:');
disp(xIntercept);

disp('Y intercept of the line:');
disp(yIntercept);

% ==============================================================================

% Part d

% Getting the connecting line between the transformed points

l2 = cross(x3, y3);
disp('Connecting line between the two points x3 and y3:'), disp(l2);

xIntercept2 = (-l2(3)/l2(1));
yIntercept2 = (-l2(3)/l2(2));

disp('X intercept of the new line:');
disp(xIntercept2);

disp('Y intercept of the new line:');
disp(yIntercept2);

% ==============================================================================

% Part e

plot([x(1) y(1)], [x(2) y(2)], 'go-', [x3(1) y3(1)], [x3(2) y3(2)], 'b*-');
xlabel('x');
ylabel('y');
