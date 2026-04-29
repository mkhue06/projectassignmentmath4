clear; clc;

A = 2;
B = -2;
C = 1;

x = 6;
y = 10;

tol = 1e-6;
max_iter = 100;

H = [2*A, -B; 
     -B, 2*C];

results = [];

for n = 1:max_iter

    fval = A*x^2 - B*x*y + C*y^2 + x - y;

    c = [2*A*x - B*y + 1;
         -B*x + 2*C*y - 1];

    grad_norm = norm(c);

    alpha = (c' * c) / (c' * H * c);

    x_new = x - alpha*c(1);
    y_new = y - alpha*c(2);

    step_size = norm([x_new - x, y_new - y]);

    results = [results; n, x, y, fval, c(1), c(2), alpha, grad_norm, step_size];

    if step_size < tol
        x = x_new;
        y = y_new;
        break;
    end

    x = x_new;
    y = y_new;

end

f_min = A*x^2 - B*x*y + C*y^2 + x - y;

T = array2table(results, ...
    'VariableNames', {'Iteration','x','y','f_value','grad_x','grad_y','alpha','grad_norm','step_size'});

disp(T)

fprintf('\nAnalytical gradient result:\n');
fprintf('x = %.10f\n', x);
fprintf('y = %.10f\n', y);
fprintf('f_min = %.10f\n', f_min);
fprintf('Iterations = %d\n', n);
