close all
clear all;
clc
a = [0 3 10];
b = [1 5 7 13];

sums = zeros(1, numel(a) * numel(b));
for i = 1:numel(b)
    for j = 1:numel(a)
        sums((i-1)*numel(a) + j) = a(j) + b(i); 
    end
end

if numel(unique(sums)) ~= numel(sums)
    error('bad choice of a,b');
end

sums= sort(sums);

diffs = zeros(1,numel(sums) * numel(sums));
for i = numel(sums):-1:1
    for j = numel(sums):-1:1
        diffs((i-1)*(numel(sums)) + j) = sums(i) - sums(j);
    end
end
diffs(diffs==0) = [];
uniq = unique(diffs);
numbers = min(uniq): 1: max(uniq)+1;
numbers = numbers - 0.5;

histogram(diffs, numbers)
hold on;
plot(a, 1, 'r*');
plot(b, 1, 'b*');
