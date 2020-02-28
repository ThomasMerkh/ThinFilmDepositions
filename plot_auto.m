auto = csvread('AutoCorrelated50000.csv');
auto = auto';

hr = log(auto(:,3));
r = log(auto(:,1));
x = linspace(0,4,512);
y = 0.52.*x - 2.67;

hold on;
plot(r,hr);
plot(x,y);

