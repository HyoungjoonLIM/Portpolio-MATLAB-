%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 2018311004 Hyoung-joon Lim %%%%%%%%%%%
%%%%%%%%%%%%%%% In my lab-computer %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.
% 1) using Fixed model
y = [-113 0.003 0.001 112.993 112.998 0.003 0.001 -112.999]';
A = [1 0.764 5.96 0 0 0; 0 0 0 1 0.764 5.96;
    1 5.062 10.541 0 0 0; 0 0 0 1 5.062 10.541; 
    1 9.663 6.243 0 0 0; 0 0 0 1 9.663 6.243; 
    1 5.35 1.654 0 0 0; 0 0 0 1 5.35 1.654];
std = [0.104 0.112 0.096 0.12 0.112 0.088 0.096 0.104];     
Q = diag(std.^2);
P = inv(Q);
psi_1_1 = inv(A'*P*A)*A'*P*y
e = y-A*psi_1_1;
var_1_1 = (e'*P*e)/(8-rank(A))
D_1_1 = var_1_1*inv(A'*P*A)

% 2) using General LESS
psi_1_2 = psi_1_1;
B_bar = [psi_1_2(2) psi_1_2(3) -1 0 0 0 0 0 0 0 0 0 0 0 0 0
    psi_1_2(5) psi_1_2(6) 0 -1 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 psi_1_2(2) psi_1_2(3) -1 0 0 0 0 0 0 0 0 0
    0 0 0 0 psi_1_2(5) psi_1_2(6) 0 -1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 psi_1_2(2) psi_1_2(3) -1 0 0 0 0 0
    0 0 0 0 0 0 0 0 psi_1_2(5) psi_1_2(6) 0 -1 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 psi_1_2(2) psi_1_2(3) -1 0
    0 0 0 0 0 0 0 0 0 0 0 0 psi_1_2(5) psi_1_2(6) 0 -1];
A_bar = [0.764 5.96 1 0 0 0; 0 0 0 0.764 5.96 1;
    5.062 10.541 1 0 0 0; 0 0 0 5.062 10.541 1;
    9.663 6.243 1 0 0 0; 0 0 0 9.663 6.243 1;
    5.35 1.654 1 0 0 0; 0 0 0 5.35 1.654 1];
std_p = [sqrt((psi_1_2(2)*std(1))^2+(psi_1_2(3)*std(2))^2)
    sqrt((psi_1_2(5)*std(1))^2+(psi_1_2(6)*std(2))^2)
    sqrt((psi_1_2(2)*std(3))^2+(psi_1_2(3)*std(4))^2)
    sqrt((psi_1_2(5)*std(3))^2+(psi_1_2(6)*std(4))^2)
    sqrt((psi_1_2(2)*std(5))^2+(psi_1_2(3)*std(6))^2)
    sqrt((psi_1_2(5)*std(5))^2+(psi_1_2(6)*std(6))^2)
    sqrt((psi_1_2(2)*std(7))^2+(psi_1_2(3)*std(8))^2)
    sqrt((psi_1_2(5)*std(7))^2+(psi_1_2(6)*std(8))^2)];
Q = diag([0.104 0.112 std_p(1) std_p(2) 0.096 0.12 std_p(3) std_p(4) 0.112 0.088 std_p(5) std_p(6) 0.096 0.104 std_p(7) std_p(8)].^2);
P1 = inv(Q);

delta = inv(A_bar'*inv(B_bar*inv(P1)*B_bar')*A_bar)*A_bar'*inv(B_bar*inv(P1)*B_bar')*e
psi_1_2 = psi_1_2 - delta;
e = y-A*psi_1_2;
var_1_2 = (e'*P*e)/(8-rank(A))
D_1_2 = var_1_2*inv(A_bar'*inv(B_bar*inv(P1)*B_bar')*A_bar)

% 3) using Error in variables
Q_x = diag(std.^2);
Q_0 = diag([1 1 0 1 1 0]);
Q_y = diag(std_p.^2);

prev = psi_1_1;
psi_1_3 = inv(A'*inv(Q_y + prev'*Q_0*prev*Q_x)*A)*A'*inv(Q_y + prev'*Q_0*prev*Q_x)*y
n=1;

while max(abs(psi_1_3-prev)) > 10^-5
    prev = psi_1_3;
    psi_1_3 = inv(A'*inv(Q_y + prev'*Q_0*prev*Q_x)*A)*A'*inv(Q_y + prev'*Q_0*prev*Q_x)*y
    n=n+1
end

e = y-A*psi_1_3;
var_1_3 = (e'*P*e)/(8-rank(A))
D_1_3 = var_1_3*inv(A'*inv(Q_y + inv(prev'*Q_0*prev*Q_x))*A)
    
% 4) 
P_1_1 = [1 1.746 9.354 0 0 0; 0 0 0 1 1.746 9.354]*psi_1_1
P_1_2 = [1 1.746 9.354 0 0 0; 0 0 0 1 1.746 9.354]*psi_1_2
P_1_3 = [1 1.746 9.354 0 0 0; 0 0 0 1 1.746 9.354]*psi_1_3

P_2_1 = [1 5.329 9.463 0 0 0; 0 0 0 1 5.329 9.463]*psi_1_1
P_2_2 = [1 5.329 9.463 0 0 0; 0 0 0 1 5.329 9.463]*psi_1_2
P_2_3 = [1 5.329 9.463 0 0 0; 0 0 0 1 5.329 9.463]*psi_1_3



%% 2.
% 1) Condition Eqns
P = diag([1,2,2]);
B = [1 1 1];
y = [14+25/60+20/3600 58+16/60 107+19/60+10/3600]';
k = 180;
ans_2_1 = y-inv(P)*B'*inv(B*inv(P)*B')*(B*y-k)

% 2) Obs Eqns
A = [1 0; 0 1; -1 -1];
y = [14+25/60+20/3600 58+16/60 107+19/60+10/3600-180]';
ans_2_2 = inv(A'*P*A)*A'*P*y;
ans_2_2 = [ans_2_2; 180-sum(ans_2_2)]

%% 3.
% 1) using Fixed model
y = [1420.407 895.362 895.887 351.398 -944.926 641.434 968.084 -1384.138 1993.262 -2367.511 -3382.284 3487.762]';
A = [90 90 1 0 0 0 -y(1)*90 -y(1)*90
     0 0 0 90 90 1 -y(2)*90 -y(2)*90
     50 40 1 0 0 0 -y(3)*50 -y(3)*40
     0 0 0 50 40 1 -y(4)*50 -y(4)*40
     -30 20 1 0 0 0 -y(5)*(-30) -y(5)*20
     0 0 0 -30 20 1 -y(6)*(-30) -y(6)*20
     50 -40 1 0 0 0 -y(7)*50 -y(7)*(-40)
     0 0 0 50 -40 1 -y(8)*50 -y(8)*(-40)
     110 -80 1 0 0 0 -y(9)*110 -y(9)*(-80)
     0 0 0 110 -80 1 -y(10)*110 -y(10)*(-80)
     -100 80 1 0 0 0 -y(11)*(-100) -y(11)*80
     0 0 0 -100 80 1 -y(12)*(-100) -y(12)*80];
P = eye(12);
psi = inv(A'*P*A)*A'*P*y
e = y-A*psi;
var_3_1 = (e'*P*e)/(12-8)
D_3_1 = var_3_1*inv(A'*P*A)

% 2) using General LESS
B_bar = []; 
A_bar = [];
xy = [90 90; 50 40; -30 20; 50 -40; 110 -80; -100 80];
XY = [1420.407 895.362; 895.887 351.398; -944.926 641.434; 
         968.084 -1384.138; 1993.262 -2367.511; -3382.284 3487.762];

syms a1 b1 c1 a2 b2 c2 a3 b3 x y X Y

F = (a1*x + b1*y + c1)/(a3*x + b3*y + 1)-X;
G = (a2*x + b2*y + c2)/(a3*x + b3*y + 1)-Y; 

FG = [F G];
JA = jacobian(FG, [a1, b1, c1, a2, b2, c2, a3, b3]); 
JB = jacobian(FG, [x, y]);    

for i=1:6
    A_bar = [A_bar; subs(JA, {x y}, {xy(i,1) xy(i,2)})];
    B_bar = [B_bar; subs(JB, {x y}, {xy(i,1) xy(i,2)})];           
end

A_bar = eval(subs(A_bar, {a1 b1 c1 a2 b2 c2 a3 b3}, psi'));
B_bar = eval(subs(B_bar, {a1 b1 c1 a2 b2 c2 a3 b3}, psi'));

temp = [-1 0; 0 -1; -1 0; 0 -1; -1 0; 0 -1; -1 0; 0 -1; -1 0; 0 -1; -1 0; 0 -1];
B_bar = [B_bar temp];

F2 = []; G2 = [];
for i=1:6
    F2 = [F2; subs(F, {x y X}, {xy(i,1) xy(i,2) XY(i,1)})];
    G2 = [G2; subs(G, {x y Y}, {xy(i,1) xy(i,2) XY(i,2)})];
end

F2 = eval(subs(F2, {a1 b1 c1 a2 b2 c2 a3 b3}, psi'));
G2 = eval(subs(G2, {a1 b1 c1 a2 b2 c2 a3 b3}, psi'));

w = [F2; G2];

B_bar1 = [B_bar(1, 1:4) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
          B_bar(2, 1:4) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 B_bar(3, 1:4) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 B_bar(4, 1:4) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 B_bar(5, 1:4) 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 B_bar(6, 1:4) 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0 0 B_bar(7, 1:4) 0 0 0 0 0 0 0 0 
          0 0 0 0 0 0 0 0 0 0 0 0 B_bar(8, 1:4) 0 0 0 0 0 0 0 0 
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 B_bar(9, 1:4) 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 B_bar(10, 1:4) 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 B_bar(11, 1:4)
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 B_bar(12, 1:4)];
std = 0.3;
std_p = [sqrt((B_bar(1,1)*std)^2+(B_bar(1,2)*std)^2)
         sqrt((B_bar(2,1)*std)^2+(B_bar(2,2)*std)^2)
         sqrt((B_bar(3,1)*std)^2+(B_bar(3,2)*std)^2)
         sqrt((B_bar(4,1)*std)^2+(B_bar(4,2)*std)^2)
         sqrt((B_bar(5,1)*std)^2+(B_bar(5,2)*std)^2)
         sqrt((B_bar(6,1)*std)^2+(B_bar(6,2)*std)^2)
         sqrt((B_bar(7,1)*std)^2+(B_bar(7,2)*std)^2)
         sqrt((B_bar(8,1)*std)^2+(B_bar(8,2)*std)^2)
         sqrt((B_bar(9,1)*std)^2+(B_bar(9,2)*std)^2)
         sqrt((B_bar(10,1)*std)^2+(B_bar(10,2)*std)^2)
         sqrt((B_bar(11,1)*std)^2+(B_bar(11,2)*std)^2)
         sqrt((B_bar(12,1)*std)^2+(B_bar(12,2)*std)^2)];
     
Q = diag([0.3 0.3 std_p(1) std_p(2) 0.3 0.3 std_p(3) std_p(4) 0.3 0.3 std_p(5) std_p(6) 0.3 0.3 std_p(7) std_p(8) 0.3 0.3 std_p(9) std_p(10) 0.3 0.3 std_p(11) std_p(12)].^2);
P1=inv(Q);

delta = inv(A_bar'*inv(B_bar1*inv(P1)*B_bar1')*A_bar)*A_bar'*inv(B_bar1*inv(P1)*B_bar1')*w;
psi_3_2 = psi-delta

y = [1420.407 895.362 895.887 351.398 -944.926 641.434 968.084 -1384.138 1993.262 -2367.511 -3382.284 3487.762]';
e = y-A*psi_3_2;
var_3_2 = (e'*P*e)/(12-8)
D_3_2 = var_3_2*inv(A_bar'*inv(B_bar1*inv(P1)*B_bar1')*A_bar)

% 4) 
xy1 = [-60 20]; xy2 = [-100 -100];
syms a1 b1 c1 a2 b2 c2 a3 b3 x y
X=(a1*x + b1*y + c1)/(a3*x + b3*y + 1);
Y=(a2*x + b2*y + c2)/(a3*x + b3*y + 1);

XY = [X Y];

P_3_1 = eval(subs(XY, {a1 b1 c1 a2 b2 c2 a3 b3}, psi_3_2'));
P_3_1 = eval(subs(P_3_1, {x y}, xy1))

P_3_2 = eval(subs(XY, {a1 b1 c1 a2 b2 c2 a3 b3}, psi_3_2'));
P_3_2 = eval(subs(P_3_2, {x y}, xy2))


%% 4. 
% 1) fixed constraint
k_0 = [0 0 0 0]';
ini = [3048 3048 3300 3400 3806.062 3209.133 3600 2800 3300 2800 3500 3100];
syms xa ya xb yb xc yc xd yd xe ye xf yf
dab = sqrt((xa-xb)^2 + (ya-yb)^2);
dae = sqrt((xa-xe)^2 + (ya-ye)^2);
dbe = sqrt((xb-xe)^2 + (yb-ye)^2);
dbf = sqrt((xb-xf)^2 + (yb-yf)^2);
dbc = sqrt((xb-xc)^2 + (yb-yc)^2);

dcf = sqrt((xc-xf)^2 + (yc-yf)^2);
dcd = sqrt((xc-xd)^2 + (yc-yd)^2);
ddf = sqrt((xd-xf)^2 + (yd-yf)^2);
dde = sqrt((xd-xe)^2 + (yd-ye)^2);
def = sqrt((xe-xf)^2 + (ye-yf)^2);

delta = 999;
n=0;

while max(abs(delta)) > 10^-4
    dist = [dab dae dbe dbf dbc dcf dcd ddf dde def];
    dist1 = eval(subs(dist, {xa ya xb yb xc yc xd yd xe ye xf yf}, ini));

    y_0 = [426.997 332.4 501.18 371.106 525.308 297.564 379.293 256.87 318.513 283.747]';
    y = dist1'- y_0;
    
    A = jacobian(dist, [xa ya xb yb xc yc xd yd xe ye xf yf]);
    A = eval(subs(A, {xa ya xb yb xc yc xd yd xe ye xf yf}, ini));

    K = [1 0 0 0 0 0 0 0 0 0 0 0
         0 1 0 0 0 0 0 0 0 0 0 0
         0 0 0 0 1 0 0 0 0 0 0 0
         0 0 0 0 0 1 0 0 0 0 0 0];
    Q = diag([0.007 0.0067 0.007 0.0067 0.007 0.0067 0.007 0.0067 0.0067 0.0067].^2);
    P=inv(Q);

    N=A'*P*A; c = A'*P*y;
    delta = inv(N+K'*K)*c + inv(N+K'*K) * K' * inv(K*inv(N+K'*K)*K') * (k_0-K*inv(N+K'*K)*c);
    n=n+1;
    ini = ini-delta';
end 
psi_4_1 = ini

% 2) stochastic constraint
ini2 = [3048 3048 3300 3400 3806.062 3209.133 3600 2800 3300 2800 3500 3100];
delta = 999;
n=0;

while max(abs(delta)) > 10^-4
    dist = [dab dae dbe dbf dbc dcf dcd ddf dde def];
    dist1 = eval(subs(dist, {xa ya xb yb xc yc xd yd xe ye xf yf}, ini2));

    y_0 = [426.997 332.4 501.18 371.106 525.308 297.564 379.293 256.87 318.513 283.747]';
    y = dist1'- y_0;
    
    A = jacobian(dist, [xa ya xb yb xc yc xd yd xe ye xf yf]);
    A = eval(subs(A, {xa ya xb yb xc yc xd yd xe ye xf yf}, ini2));

    K = [1 0 0 0 0 0 0 0 0 0 0 0
         0 1 0 0 0 0 0 0 0 0 0 0
         0 0 0 0 1 0 0 0 0 0 0 0
         0 0 0 0 0 1 0 0 0 0 0 0];
    Q = diag([0.007 0.0067 0.007 0.0067 0.007 0.0067 0.007 0.0067 0.0067 0.0067].^2);
    P=inv(Q);
    Q_0 = diag([0.003 0.003 0.003 0.003].^2);
    P_0 = inv(Q_0);

    N=A'*P*A; c = A'*P*y;
    delta = inv(N+K'*K)*c + inv(N+K'*K) * K' * inv(P_0 + K*inv(N+K'*K)*K') * (k_0 - K*inv(N+K'*K)*c);
    n=n+1;
    ini2 = ini2-delta';
end 
psi_4_2 = ini2

%% 5. 
% 1) fixed constraint
k_0 = [0];
ini = [1003.07 3640 2323.07 3638.46 2496.08 1061.74];

syms xa ya xb yb xc yc xd yd
dac = sqrt((1000-xc)^2 + (1000-yc)^2);
dab = sqrt((1000-xb)^2 + (1000-yb)^2);

dbc = sqrt((xb-xc)^2 + (yb-yc)^2);
dcd = sqrt((xc-xd)^2 + (yc-yd)^2);
dda = sqrt((xd-1000)^2 + (yd-1000)^2);
dbd = sqrt((xb-xd)^2 + (yb-yd)^2);

dir = atand((xb-1000)/(yb-1000));

delta = 999;
n=0;

while max(abs(delta)) > 10^-4
    dist = [dac dab dbc dcd dda dbd];
    dist1 = eval(subs(dist, {xb yb xc yc xd yd}, ini));
    
    y_0 = [2951.604 2640.017 1320.016 2582.534 1497.360 2979.325]';
    y = dist1'- y_0;
    
    A = jacobian(dist, [xb yb xc yc xd yd]);
    A = eval(subs(A, {xb yb xc yc xd yd}, ini));

    K = jacobian(dir, [xb yb xc yc xd yd]);
    K = eval(subs(K, {xb yb xc yc xd yd}, ini));
    Q = diag([0.025 0.024 0.021 0.024 0.021 0.025].^2);
    P=inv(Q);

    N=A'*P*A; c = A'*P*y;
    delta = inv(N+K'*K)*c + inv(N+K'*K) * K' * inv(K*inv(N+K'*K)*K') * (k_0-K*inv(N+K'*K)*c);
    n=n+1;
    ini = ini-delta';
end 
psi_5_1 = ini
dir = atand((psi_5_1(1)-1000)/(psi_5_1(2)-1000))

% 2) 
clear;
ini2 = [1003.07 3640 2323.07 3638.46 2496.08 1061.74];
k_0 = [0];
syms xa ya xb yb xc yc xd yd
dac = sqrt((1000-xc)^2 + (1000-yc)^2);
dab = sqrt((1000-xb)^2 + (1000-yb)^2);

dbc = sqrt((xb-xc)^2 + (yb-yc)^2);
dcd = sqrt((xc-xd)^2 + (yc-yd)^2);
dda = sqrt((xd-1000)^2 + (yd-1000)^2);
dbd = sqrt((xb-xd)^2 + (yb-yd)^2);

dir = atand((xb-1000)/(yb-1000));

delta = 999;
n=0;
while max(abs(delta)) > 10^-4
    dist = [dac dab dbc dcd dda dbd];
    dist1 = eval(subs(dist, {xb yb xc yc xd yd}, ini2));
    
    y_0 = [2951.604 2640.017 1320.016 2582.534 1497.360 2979.325]';
    y = dist1'- y_0;
    
    A = jacobian(dist, [xb yb xc yc xd yd]);
    A = eval(subs(A, {xb yb xc yc xd yd}, ini2));

    K = jacobian(dir, [xb yb xc yc xd yd]);
    K = eval(subs(K, {xb yb xc yc xd yd}, ini2));
    Q = diag([0.025 0.024 0.021 0.024 0.021 0.025].^2);
    P=inv(Q);
    Q_0 = 10/3600^2;
    P_0 = inv(Q_0);

    N=A'*P*A; c = A'*P*y;
    delta = inv(N+K'*K)*c + inv(N+K'*K) * K' * inv(P_0 + K*inv(N+K'*K)*K') * (k_0 - K*inv(N+K'*K)*c);
    n=n+1;
    ini2 = ini2-delta';
end 
psi_5_2 = ini2
dir = atand((psi_5_2(1)-1000)/(psi_5_2(2)-1000))

















