W=2; 
H=1;
e = 0.1;
d = 0.02;


Point(1)={0,0,0};
Point(2)={W,0,0};
Point(3)={W,H,0};
Point(4)={0,H,0};
Point(5)={-e,H/2-e/2,0};
Point(6)={-e,H/2+e/2,0};
Point(7)={0,H/2-e/2,0};
Point(8)={0,H/2+e/2,0};
Point(9)={W+e,H/2-e/2,0};
Point(10)={W+e,H/2+e/2,0};
Point(11)={W,H/2-e/2,0};
Point(12)={W,H/2+e/2,0};
Point(13)={W/2-d/2,0,0};
Point(14)={W/2+d/2,0,0};
Point(15)={W/2-d/2,H/2+e,0};
Point(16)={W/2+d/2,H/2+e,0};

Point(17)={W/4-d/2,H,0};
Point(18)={W/4+d/2,H,0};
Point(19)={W/4-d/2,H/2-e,0};
Point(20)={W/4+d/2,H/2-e,0};

Point(21)={3*W/4-d/2,H,0};
Point(22)={3*W/4+d/2,H,0};
Point(23)={3*W/4-d/2,H/2-e,0};
Point(24)={3*W/4+d/2,H/2-e,0};


Line(1) = {5, 6};
Line(2) = {8, 6};
Line(3) = {7, 5};
Line(4) = {1, 7};
Line(5) = {8, 4};
Line(6) = {10, 9};
Line(7) = {9, 11};
Line(8) = {11, 2};
Line(9) = {10, 12};
Line(10) = {12, 3};
Line(11) = {4, 17};
Line(12) = {17, 19};
Line(13) = {19, 20};
Line(14) = {20, 18};
Line(15) = {18, 21};
Line(16) = {21, 23};
Line(17) = {23, 24};
Line(18) = {24, 22};
Line(19) = {22, 3};
Line(20) = {2, 14};
Line(21) = {14, 16};
Line(22) = {16, 15};
Line(23) = {15, 13};
Line(24) = {13, 1};

Line Loop(25) = {1, -2, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, -10, -9, 6, 7, 8, 20, 21, 22, 23, 24, 4, 3};

Plane Surface(26) = {25};
Physical Surface(27) = {26};

//INLET 
Physical Line(71) = {1};
//OUTLET 
Physical Line(72) = {6};

//WALLS
Physical Curve(73) = {2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};





Mesh.CharacteristicLengthFactor=0.3;
Mesh.Algorithm=1;
Mesh.Format=30;
Mesh.ScalingFactor = 1.0;
Mesh.Smoothing = 10;
