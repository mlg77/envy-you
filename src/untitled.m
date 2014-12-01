% %  plot getRef input:
startpulse = 20;
lengthpulse = 50;

for i = 1:100
    tx(i) = i;
    y(i) = getRef2(i,'ft');
end
figure; plot(t,y)

for i = 1:100
    tx(i) = i;
    y(i) = getRef2(i,'fluxft');
end
figure; plot(t,y)