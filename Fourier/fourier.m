clf;
clc;
clear;
lambda = 0.24;
A = 0.44;
P = 1.0;
W = 0.5;
n = 100;
k = 20;
X = linspace(-2*P,2*P);
a0 = 0.5;
b1 = integral(@(x) 2/P*sin(2*pi*x/P), 0, P/2);
% 
f = cell(k,1);
f{1,1} = @(x) (a0 + b1.*sin(2*pi*x/P));

for i = 2:k
    bi = integral(@(x) 2/P*sin(2*pi*i*x/P), 0, P/2);
    if (bi ~= 0)
      f{i,1} = @(x) (bi*sin(2*pi*i*x/P));
    else
        continue;
    end
    
end
figure(1);
hold on
Y = sin(X*2*pi/P);
Y(find(Y >= 0)) = 1; 
Y(find(Y < 0)) = 0;
fsum = @(x)0;
for i = 1:k
    clf;
    hold on;
    Y = sin(X*2*pi/P);
    Y(find(Y >= 0)) = 1; 
    Y(find(Y < 0)) = 0;
    plot(X,Y);
    fsum = @(x) fsum(x) + f{i}(x);
    plot (X,fsum(X));
    pause(0.5);
end
title('График распределения интенсивности');
xlabel('Координата');
ylabel('Интенсивность');
nlim = floor(2*A/lambda);
counter = 0;
while (counter*lambda / (2 * A) <= 1)
    R = counter*lambda / (2 * A);   
    counter = counter + 1;
end
stepper = 0:(counter - 1);
R = stepper*lambda / (2 * A); 
T = 2 / pi * (acos(R) - R.*sqrt(1-R.^2));
figure(2);
plot(stepper,T);
title('Функция передачи модуляции');
xlabel('Период на единицу длины');
ylabel('КМП');