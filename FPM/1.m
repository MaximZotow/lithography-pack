%clf;

clc;
clear;
lambda = 0.24;
A = 0.44;
P = 1;
W = P / 2;
n = 1000;
nlim = floor(2*A*P/lambda);
nki=1:nlim;

R = nki / P *lambda / (2 * A); 

T = 2 / pi * (acos(R) - R.*sqrt(1-R.^2));
%FIG1
figure(1);
%FIG1
clf;
hold on;
grid on;
xlabel('x');
ylabel('I_0');
axis([-P/2 P/2 -0.2 1.2]);
X = -P / 2:P/100:P / 2;
    Y = sin(X*2*pi/P);
    Y(find(Y >= 0)) = 1; 
    Y(find(Y < 0)) = 0;
    axis([-P/2 P/2 -0.2 1.2]);
plot(X,Y);
axis([-P/2 P/2 -0.2 1.2]);
filename_fig1 = "fig1-0.pdf";
  %print (figure(1), filename_fig1);
  

a0 = 0.5;
%b1 = integral(@(x) 2/P*sin(2*pi*x/P), 0, P/2);

f = cell(nki,1);
for i=nki
  bi = integral(@(x) 2/P*sin(2*pi*i*x/P), 0, P/2);
  f{i,1}=@(x) bi.*sin(2*pi*x*i/P);
end
%FIG1
figure(1);
%FIG1
xlabel('x');
ylabel('I_0');

sum_func = @(x) 0;
for i = nki
  hold on;
  grid on;
  buf = @(x) f{i}(x);
  sum_func = @(x) sum_func(x) + buf(x);
  plot(X, (buf(X) + a0));
  filename_fig = strcat("fig1-", int2str(i));
  filename_fig1 = strcat(filename_fig,".pdf");
  %print (figure(1), filename_fig1);
  pause(0.5);
end

plot(X, sum_func(X)+a0);

filename_fig = strcat("fig1-", int2str(max(nki)+1));
filename_fig1 = strcat(filename_fig, ".pdf");
%print (figure(1), filename_fig1);
fsum = @(x) a0;
%FIG2
figure(2);
%FIG2

clf;
for i = nki
    hold on;
    grid on;
    xlabel('x');
    ylabel( 'I_0');
    Y = sin(X*2*pi/P);
    Y(find(Y >= 0)) = 1; 
    Y(find(Y < 0)) = 0;
    axis([-P/2 P/2 0 1.2]);
    plot(X,Y);
    fsum = @(x) fsum(x) + T(i)*f{i}(x);
    axis([-P/2 P/2 0 1.2]);
    plot (X,fsum(X));
    filename_fig = strcat("fig2-",int2str(i));
    filename_fig2 = strcat(filename_fig, ".pdf");
    %print (figure(2), filename_fig2);
    pause(0.5);
end

fprintf("Max contrast = %d\n", max(fsum(X)));
fprintf("Min contrast = %d\n", min(fsum(X)));
contrast_ratio = (max(fsum(X)) - min(fsum(X))) / (max(fsum(X)) + min(fsum(X)));
fprintf("Contrast ratio = %d\n", contrast_ratio);


%FIG3
figure(3);
%FIG3
clf;
nki1=0:nlim;
R1 = nki1 / P *lambda / (2 * A); 
T1 = 2 / pi * (acos(R1) - R1.*sqrt(1-R1.^2));
drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),-1 );
for i=nki
  hold on;
  grid on;
  iks=[R(i) R(i)];
  ygrek=[0 T(i)];
  drawArrow(iks, ygrek);
  iks1=[R(i) 0];
  ygrek1=[T(i) T(i)];
  drawArrow(iks1, ygrek1);
  fprintf("R=%5.3f, T=%5.3f\n", R(i), T(i));
  A=max(T(i)*f{i}(X))-min(T(i)*f{i}(X));
  fprintf("Amplitude for %d harmonic A=%5.3f\n", i, A);
end
hold on;
grid on;
plot(R1,T1);
grid on;
xlabel('R');
ylabel('T');
axis([0 1 0 1]);
filename_fig3 = "func_mod.pdf";
print (figure(3), filename_fig3);
%
figure(4);
clf;
hold on;
for i=nki
  xlabel('x');
    ylabel( 'I_0');
  plot(X,f{i,1}(X)*T(i))
end;
xlabel('x');
    ylabel( 'I_0');
filename_fig4 = "fig4-1.pdf";
print (figure(4), filename_fig4);


figure(5);
clf;
res = @(x) 0.5;
for i=nki
  res = @(x) res(x) + f{i,1}(x);
end
hold on;
grid on;
xlabel('x');
    ylabel( 'I_0');
plot(X,Y);
plot(X,res(X));
filename_fig5 = "fig5-1.pdf";
print (figure(5), filename_fig5);