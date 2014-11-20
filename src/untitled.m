% %  plot getRef input:
startpulse = 20;
lengthpulse = 50;

for i = 1:100
    t(i) = i;
    y(i) = getRef(i,'ft');
end
figure; plot(t,y)

for i = 1:100
    t(i) = i;
    y(i) = getRef(i,'fluxft');
end
figure; plot(t,y)