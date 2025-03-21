clear;
clc;

% 输入D1, D2, L3', f'的值
D1 = -1200;
D2 = 1200;
L3_prime = -1300;
f_prime =-2000;

% 定义符号变量
syms alpha1 alpha2 beta1 beta2

% 定义四个方程
eqn1 = D1 - (1 - alpha1)/(beta1*beta2)*f_prime == 0;
eqn2 = D2 - alpha1*(1 - alpha2)/beta2*f_prime == 0;
eqn3 = L3_prime - alpha1*alpha2*f_prime == 0;
eqn4 = beta1*beta2 - beta2*(1 + beta1)/alpha1 + (1 + beta2)/(alpha1*alpha2) == 0;

% 求解方程组
sol = solve([eqn1, eqn2, eqn3, eqn4], [alpha1, alpha2, beta1, beta2]);

% 尝试将解转换为数值形式
try
    alpha1_sol = double(sol.alpha1);
    alpha2_sol = double(sol.alpha2);
    beta1_sol = double(sol.beta1);
    beta2_sol = double(sol.beta2);
    disp('alpha1的数值解:');
    disp(alpha1_sol);
    disp('alpha2的数值解:');
    disp(alpha2_sol);
    disp('beta1的数值解:');
    disp(beta1_sol);
    disp('beta2的数值解:');
    disp(beta2_sol);
catch
    disp('无法将解转换为数值形式，可能是没有实数解或者解过于复杂。');
    disp('符号解如下:');
    disp('alpha1的解:');
    disp(sol.alpha1);
    disp('alpha2的解:');
    disp(sol.alpha2);
    disp('beta1的解:');
    disp(sol.beta1);
    disp('beta2的解:');
    disp(sol.beta2);
end
   

a1=alpha1_sol(2);
%这有两个解可选1 2
a2=L3_prime/(a1*f_prime);
b2=(a1*(1 - a2)/D2)*f_prime;
b1=((1-a1)/(D1*b2))*f_prime;

F_prime=f_prime;
   % 计算曲率半径
    R1 = (2 / (b1 * b2)) * F_prime;
    R2 = (2 * a1 / (b2 * (1 + b1))) * F_prime;
    R3 = (2 * a1 * a2 / (1 + b2)) * F_prime;

    % 输出结果
    fprintf('遮拦比:\n a1 = %.4f\n a2 = %.4f\n', a1, a2);
    fprintf('放大率:\n b1 = %.4f\n b2 = %.4f\n', b1, b2);
    fprintf('曲率半径（mm）:\n R1 = %.4f\n R2 = %.4f\n R3 = %.4f\n', R1, R2, R3);

    
    



% 定义 s1 函数
function s1 = calculate_s1(e1, e2, e3, a1, a2, b1, b2)
    term1 = 0.25 * ((e1 - 1) * b1^3 * b2^3 - e2 * a1 * b2^3 * (1 + b1)^3);
    term2 = e3 * a1 * a2 * (1 + b2)^3;
    term3 = a1 * b2^3 * (1 + b1) * (1 - b1)^2;
    term4 = -a1 * a2 * (1 + b2) * (1 - b2)^2;
    s1 = term1 + term2 + term3 + term4;
end

% 定义 s2_prime 函数
function s2_prime = calculate_s2_prime(e1, e3, a1, a2, b1, b2)
    term1 = ((1 - a1) * b1^2 * b2^2 * e1) / (4 * a1);
    term2 = ((1 - a1) * (1 + b2)^3 * e3) / (4 * b1);
    term3 = -((1 - a1) * b1^2 * b2^2) / (4 * a1);
    term4 = -((1 - a2) * (1 + b1) * (1 - b2)^2) / (4 * b2);
    term5 = -0.5;
    s2_prime = term1 + term2 + term3 + term4 + term5;
end

% 定义 s3_prime 函数
function s3_prime = calculate_s3_prime(e1, e3, a1, a2, b1, b2)
    term1 = e1 * ((a1 - 1)^2 * b1 * b2) / (4 * a1^2);
    term2 = e3 * ((1 - a2)^2 * (1 + b2)^3) / (4 * a1 * a2 * b2^2);
    term3 = -((a1 - 1)^2 * b1 * b2) / (4 * a1^2);
    term4 = -((1 - a2)^2 * (1 + b2) * (1 - b2)^2) / (4 * a1 * a2 * b2^2);
    term5 = -((1 - a1) * b1 * b2) / a1;
    term6 = -((1 - a2) * (1 - b2)^2) / (a1 * a2 * b2);
    term7 = -b1 * b2;
    term8 = b2 * (1 + b1) / a1;
    term9 = -(1 + b2) / (a1 * a2);
    s3_prime = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9;
end

% 定义方程组函数
function F = equations(x, a1, a2, b1, b2)
    e1 = x(1);
    e2 = x(2);
    e3 = x(3);
    F(1) = calculate_s1(e1, e2, e3, a1, a2, b1, b2);
    F(2) = calculate_s2_prime(e1, e3, a1, a2, b1, b2);
    F(3) = calculate_s3_prime(e1, e3, a1, a2, b1, b2);
end


% 初始猜测值
x0 = [1; 1; 1];

% 使用 fsolve 求解方程组
[x, fval] = fsolve(@(x) equations(x, a1, a2, b1, b2), x0);

% 提取解
e1 = x(1);
e2 = x(2);
e3 = x(3);

% 输出结果
fprintf('求解得到的 e1 = %.6f, e2 = %.6f, e3 = %.6f\n', e1, e2, e3);
fprintf('方程组的残差: fval = [%.6f, %.6f, %.6f]\n', fval(1), fval(2), fval(3));    
