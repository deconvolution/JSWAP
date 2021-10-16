tt=readtable('receiver_list.txt');
range_lon=[16.0761, 16.4618];
range_lat=[39.1452,39.5106];
range_z=[-600,18000];
for i=1:24
    if 1-ismember(i,[1,10,17])
    tt=load(['./event_data/event_data_',num2str(i),'.mat']);
    
    if str2double(tt.data.S{4})<range_lat(1)
        range_lat(1)=str2double(tt.data.S{4});
    end
    if str2double(tt.data.S{4})>range_lat(2)
        range_lat(2)=str2double(tt.data.S{4});
    end
    if str2double(tt.data.S{5})<range_lon(1)
        range_lon(1)=str2double(tt.data.S{5});
    end
    if str2double(tt.data.S{5})>range_lon(2)
        range_lon(2)=str2double(tt.data.S{5});
    end
    if str2double(tt.data.S{6})*1000>range_z(2)
        range_z(2)=str2double(tt.data.S{6})*1000;
        disp(i)
    end
    end
end
range_lon
range_lat
range_z
%%
dx=150;
dy=150;
dz=150;
%% for input
listing=dir('./event_data');
for i=1:size(listing,1)-2
    %% generate empty file
    data.Rp=[];
    data.Rs=[];
    data.S=[];
    tt=load(listing(i+2).name);
    %% location
    data.Rp=fix((cell2mat(tt.data.Rp(:,3))-range_lon(1))*86000/dx);
    data.Rp=[data.Rp,fix((cell2mat(tt.data.Rp(:,2))-range_lat(1))*111000/dy)];
    data.Rp=[data.Rp,fix((cell2mat(tt.data.Rp(:,4))-range_z(1))/dz)];
    
    data.Rs=fix((cell2mat(tt.data.Rs(:,3))-range_lon(1))*86000/dx);
    data.Rs=[data.Rs,fix((cell2mat(tt.data.Rs(:,2))-range_lat(1))*111000/dy)];
    data.Rs=[data.Rs,fix((cell2mat(tt.data.Rs(:,4))-range_z(1))/dz)];
    
    data.S=fix((str2num(cell2mat(tt.data.S(5)))-range_lon(1))*86000/dx);
    data.S=[data.S,fix((str2num(cell2mat(tt.data.S(4)))-range_lat(1))*111000/dy)];
    data.S=[data.S,fix((str2num(cell2mat(tt.data.S(6)))*1000-range_z(1))/dz)];
    %% time
    data.Rp=[data.Rp,cell2mat(tt.data.Rp(:,6))-str2num(cell2mat(tt.data.S(2)))];
    data.Rs=[data.Rs,cell2mat(tt.data.Rs(:,6))-str2num(cell2mat(tt.data.S(2)))];
    for j=1:size(data.Rp,1)
        data.Rp(j,4)=data.Rp(j,4)+seconds(tt.data.Rp{j,7}-duration(tt.data.S{:,3}));
        if data.Rp(j,4)<0
            break
            disp('error');
        end
    end
    for j=1:size(data.Rs,1)
         data.Rs(j,4)=data.Rs(j,4)+seconds(tt.data.Rs{j,7}-duration(tt.data.S{:,3}));
         if data.Rs(j,4)<0
             break
             disp('error');
         end
    end
    save(['./crati_traveltime_input/' listing(i+2).name],'data');
end