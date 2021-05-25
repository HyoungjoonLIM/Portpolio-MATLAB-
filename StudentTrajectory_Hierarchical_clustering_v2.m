%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Implementation of Hierarchical Clustering %%%%
%%%%%%%%%%%%% 20180422 Hyoung-joon Lim %%%%%%%%%%%%
%%%%%%%%%%%%%%% In my lab-computer %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

cd('E:\학생데이터 분석\2017\DTW_SEM(180418)\result_dtw\201905_all')

% File reading and preprocessing for HC
s = fileread('distmat161_all.csv');
s = str2num(s);
s = str2double(s);
d_origin = [];
d_origin = s(1,:);
d_origin = [d_origin; zeros(1,size(s,2))];
d_origin = [d_origin; s(2:size(s,2),:)];
d_origin(2,1) = d_origin(1,2);
d_origin = [d_origin zeros(size(s,2)+1,1)];
d_origin(1,size(s,2)+1) = d_origin(size(s,2)+1,1);

% d_origin = distmat161natfem1234;
d1 = d_origin(2:size(d_origin,1), 2:size(d_origin,2));
for c = 2:size(d1,1) % construction of distance matrix in symmetric form
    c
    for r = 1:size(d1,1)
        if c>r
            d1(r,c) = d1(c,r); 
        end
    end
end
c1 = evalclusters(d1, 'linkage', 'Silhouette', 'KList', [2:10]);
c2 = evalclusters(d1, 'linkage', 'DaviesBouldin', 'KList', [2:10]);
c3 = evalclusters(d1, 'linkage', 'gap', 'KList', [2:10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = fileread('distmat162_natfem1234.csv');
s = str2num(s);
s = str2double(s);
d_origin = [];
d_origin = s(1,:);
d_origin = [d_origin; zeros(1,size(s,2))];
d_origin = [d_origin; s(2:size(s,2),:)];
d_origin(2,1) = d_origin(1,2);
d_origin = [d_origin zeros(size(s,2)+1,1)];
d_origin(1,size(s,2)+1) = d_origin(size(s,2)+1,1);

% d_origin = distmat162natfem1234;
d2 = d_origin(2:size(d_origin,1), 2:size(d_origin,2));
for c = 2:size(d2,1) % construction of distance matrix in symmetric form
    c
    for r = 1:size(d2,1)
        if c>r
            d2(r,c) = d2(c,r); 
        end
    end
end
c4 = evalclusters(d2, 'linkage', 'Silhouette', 'KList', [2:10]);
c5 = evalclusters(d2, 'linkage', 'DaviesBouldin', 'KList', [2:10]);
c6 = evalclusters(d2, 'linkage', 'gap', 'KList', [2:10]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = fileread('distmat171_natfem1234.csv');
s = str2num(s);
s = str2double(s);
d_origin = [];
d_origin = s(1,:);
d_origin = [d_origin; zeros(1,size(s,2))];
d_origin = [d_origin; s(2:size(s,2),:)];
d_origin(2,1) = d_origin(1,2);
d_origin = [d_origin zeros(size(s,2)+1,1)];
d_origin(1,size(s,2)+1) = d_origin(size(s,2)+1,1);

% d_origin = distmat171natfem1234;
d3 = d_origin(2:size(d_origin,1), 2:size(d_origin,2));
for c = 2:size(d3,1) % construction of distance matrix in symmetric form
    c
    for r = 1:size(d3,1)
        if c>r
            d3(r,c) = d3(c,r); 
        end
    end
end
c7 = evalclusters(d3, 'linkage', 'Silhouette', 'KList', [2:10]);
c8 = evalclusters(d3, 'linkage', 'DaviesBouldin', 'KList', [2:10]);
c9 = evalclusters(d3, 'linkage', 'gap', 'KList', [2:10]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = fileread('distmat172_natfem1234.csv');
s = str2num(s);
s = str2double(s);
d_origin = [];
d_origin = s(1,:);
d_origin = [d_origin; zeros(1,size(s,2))];
d_origin = [d_origin; s(2:size(s,2),:)];
d_origin(2,1) = d_origin(1,2);
d_origin = [d_origin zeros(size(s,2)+1,1)];
d_origin(1,size(s,2)+1) = d_origin(size(s,2)+1,1);

% d_origin = distmat172libmale1234;
d4 = d_origin(2:size(d_origin,1), 2:size(d_origin,2));
for c = 2:size(d4,1) % construction of distance matrix in symmetric form
    c
    for r = 1:size(d4,1)
        if c>r
            d4(r,c) = d4(c,r); 
        end
    end
end
c10 = evalclusters(d4, 'linkage', 'Silhouette', 'KList', [2:10]);
c11 = evalclusters(d4, 'linkage', 'DaviesBouldin', 'KList', [2:10]);
c12 = evalclusters(d4, 'linkage', 'gap', 'KList', [2:10]);


res = [c1.CriterionValues; c2.CriterionValues; c3.CriterionValues;
    c4.CriterionValues; c5.CriterionValues; c6.CriterionValues;
    c7.CriterionValues; c8.CriterionValues; c9.CriterionValues;
    c10.CriterionValues; c11.CriterionValues; c12.CriterionValues];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = clusteringindexavg;
a(:,10) = [];

% Silhouette
figure
plot((2:10),a(1,:),'k') % all = black
hold on
plot((2:10),a(5,:),'y') % engmale = yellow
hold on
plot((2:10),a(9,:),'m') % engfem = m
hold on
plot((2:10),a(13,:),'c') % libmale = cyan
hold on
plot((2:10),a(17,:),'r') % libfem = red
hold on
plot((2:10),a(21,:),'g') % natmale = green
hold on
plot((2:10),a(25,:),'b') % natfem = blue
xlabel('the number of clusters')
ylabel('Silhouette index')
legend('All students','Male students in engineering college','Female students in engineering college','Male students in liberal arts ','Female students in liberal arts','Male students in science and engineering','Female students in science and engineering','Location','northeast')

% Davies-Bouldin
figure
plot((2:10),a(2,:),'k') % all = black
hold on
plot((2:10),a(6,:),'y') % engmale = yellow
hold on
plot((2:10),a(10,:),'m') % engfem = m
hold on
plot((2:10),a(14,:),'c') % libmale = cyan
hold on
plot((2:10),a(18,:),'r') % libfem = red
hold on
plot((2:10),a(22,:),'g') % natmale = green
hold on
plot((2:10),a(26,:),'b') % natfem = blue
xlabel('the number of clusters')
ylabel('Davies-Bouldin index')
legend('All students','Male students in engineering college','Female students in engineering college','Male students in liberal arts ','Female students in liberal arts','Male students in science and engineering','Female students in science and engineering','Location','southeast')

% delta Gap index
figure
plot((2:10),a(4,:),'k') % all = black
hold on
plot((2:10),a(8,:),'y') % engmale = yellow
hold on
plot((2:10),a(12,:),'m') % engfem = m
hold on
plot((2:10),a(16,:),'c') % libmale = cyan
hold on
plot((2:10),a(20,:),'r') % libfem = red
hold on
plot((2:10),a(24,:),'g') % natmale = green
hold on
plot((2:10),a(28,:),'b') % natfem = blue
xlabel('the number of clusters')
ylabel('Gap(k+1) - Gap(k)')
legend('All students','Male students in engineering college','Female students in engineering college','Male students in liberal arts ','Female students in liberal arts','Male students in science and engineering','Female students in science and engineering','Location','northeast')


% Gap index
figure
plot((2:10),a(3,:),'k') % all = black
figure
plot((2:10),a(7,:),'y') % engmale = yellow
hold on
plot((2:10),a(11,:),'m') % engfem = m
hold on
plot((2:10),a(15,:),'c') % libmale = cyan
hold on
plot((2:10),a(19,:),'r') % libfem = red
hold on
plot((2:10),a(23,:),'g') % natmale = green
hold on
plot((2:10),a(27,:),'b') % natfem = blue
legend('All students','Male students in engineering college','Female students in engineering college','Male students in liberal arts ','Female students in liberal arts','Male students in science and engineering','Female students in science and engineering')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


z = linkage(d1,'ward'); % HC, ward linkage
cophenet(z,d1)
figure
dendrogram(z)
c_3 = cluster(z, 'maxclust', 3);
c_4 = cluster(z, 'maxclust', 4);
c_5 = cluster(z, 'maxclust', 5);
c_6 = cluster(z, 'maxclust', 6);
c_7 = cluster(z, 'maxclust', 7);
c_8 = cluster(z, 'maxclust', 8);

result = [d_origin(2:size(d_origin,1),1) c_3 c_4 c_5 c_6 c_7 c_8];
dlmwrite("clustering_result172_all(c=3_8).csv", result, 'delimiter', ',', 'precision', 15);

% hold on
% z2 = linkage(d,'complete'); % HC, ward linkage
% dendrogram(z2)

% c_2 = cluster(z, 'maxclust', 2);


% I = inconsistent(z) % 불일치계수(clustering constraint의 또 다른 종류)
% c = cluster(z, 'cutoff', 0.8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
d_origin = distmat171;
d1 = d_origin(2:size(d_origin,1), 2:size(d_origin,2));

for c = 2:size(d1,1) % construction of distance matrix in symmetric form
    c
    for r = 1:size(d1,1)
        if c>r
            d1(r,c) = d1(c,r); 
        end
    end
end
        
z = linkage(d1,'ward'); % HC, ward linkage
dendrogram(z)

c_4 = cluster(z, 'maxclust', 4);

result = [d_origin(2:size(d_origin,1),1) c_4];
csvwrite("clustering_result161(c=4).csv", result);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
d_origin = distmat172;
d1 = d_origin(2:size(d_origin,1), 2:size(d_origin,2));

for c = 2:size(d1,1) % construction of distance matrix in symmetric form
    c
    for r = 1:size(d1,1)
        if c>r
            d1(r,c) = d1(c,r); 
        end
    end
end
        
z = linkage(d1,'ward'); % HC, ward linkage
dendrogram(z)

c_4 = cluster(z, 'maxclust', 4);
c_5 = cluster(z, 'maxclust', 5);
result = [d_origin(2:size(d_origin,1),1) c_5];
csvwrite("clustering_result161(c=4).csv", result);
