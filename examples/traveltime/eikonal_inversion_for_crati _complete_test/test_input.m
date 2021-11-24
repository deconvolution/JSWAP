tt=dir('./crati_traveltime_input/');

for i=3:199
    tt2=load(['./crati_traveltime_input/' tt(i).name]);
    if tt2.data.S(3)>=124
        disp(tt(i))
    end
end