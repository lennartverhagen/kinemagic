A = time2logic(artifact{1}{1},data.time{1});
B = time2logic(artifact{2}{1},data.time{1});

figure;
plot(data.time{1},A,'b'); hold on
plot(data.time{1},B,'r');