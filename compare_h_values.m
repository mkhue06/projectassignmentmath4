clear; clc;

A = 2;
B = -2;
C = 1;

h_values = [1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10];

tol = 1e-6;
max_iter = 100;

H = [2*A, -B; 
     -B, 2*C];

comparison_results = [];

for k = 1:length(h_values)

    h = h_values(k);

    x = 6;
    y = 10;

    for n = 1:max_iter

        fval = A*x^2 - B*x*y + C*y^2 + x - y;

        fx = (A*(x+h)^2 - B*(x+h)*y + C*y^2 + (x+h) - y - fval) / h;

        fy = (A*x^2 - B*x*(y+h) + C*(y+h)^2 + x - (y+h) - fval) / h;

        c = [fx; fy];

        alpha = (c' * c) / (c' * H * c);

        x_new = x - alpha*c(1);
        y_new = y - alpha*c(2);

        step_size = norm([x_new - x, y_new - y]);

        if step_size < tol
            x = x_new;
            y = y_new;
            break;
        end

        x = x_new;
        y = y_new;

    end

    f_min = A*x^2 - B*x*y + C*y^2 + x - y;

    comparison_results = [comparison_results; h, x, y, f_min, n];

end

T = array2table(comparison_results, ...
    'VariableNames', {'h','x_star','y_star','f_min','Iterations'});

disp(T)

fprintf('\nComparison completed.\n');
