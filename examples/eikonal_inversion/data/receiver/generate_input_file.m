listing=dir('./');
%% read receiver list
receiver_list=readtable('./receiver_list.txt');
%%
I=3:26;
source_name=cell(length(I),1);
for i=I
    source_name{i-I(1)+1}=listing(i).name;
end
%%
for i=1:size(source_name,1)
     source_name{i}=extractBefore(source_name{i},".");
end
%% read source
listing2=dir('../source/');
I=3:26;
source_name2=cell(length(I),1);
for i=I
    source_name2{i-I(1)+1}=listing(i).name;
end
%%
for i=1:size(source_name2,1)
     source_name2{i}=extractBefore(source_name2{i},".");
end
%%
n=1;
for i=1:size(source_name,1)
    tt=readtable(['./',listing(i+2).name]);
    indp=tt{:,2}=="P1";
    inds=tt{:,2}=="S1";
    Rp=cell(sum(indp),7);
    Rp(:,1)=table2cell(tt(indp,1));
    Rp(:,5)=table2cell(tt(indp,3));
    Rp(:,6)=table2cell(tt(indp,4));
    Rp(:,7)=table2cell(tt(indp,5));  
   
    for j=1:sum(indp,1)
        tt2=find(contains(receiver_list{:,1},tt{j,1}));
        Rp{j,2}=receiver_list{j,2};
        Rp{j,3}=receiver_list{j,3};
        Rp{j,4}=receiver_list{j,4};
    end
    
    Rs=cell(sum(inds),7);
    Rs(:,1)=table2cell(tt(inds,1));
    Rs(:,5)=table2cell(tt(inds,3));
    Rs(:,6)=table2cell(tt(inds,4));
    Rs(:,7)=table2cell(tt(inds,5));
    for j=1:sum(inds,1)
        tt2=find(contains(receiver_list{:,1},tt{j,1}));
        Rs{j,2}=receiver_list{j,2};
        Rs{j,3}=receiver_list{j,3};
        Rs{j,4}=receiver_list{j,4};
    end
    %% find corresponding source
    tt3=find(contains(source_name2,source_name{i}));
    %%
    fid=fopen(['../source/',listing2(i+2).name]);
    cell_data=textscan(fid,'%s','Delimiter',' ','headerLines',0);
    fclose(fid);
    data.S=cell_data{1}';
    data.Rp=Rp;
    data.Rs=Rs;
    save(['./event_data/event_data_' num2str(n) '.mat'],'data');
    n=n+1;
end