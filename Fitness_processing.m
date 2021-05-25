% 2017-07-17
% Lim, Hyoung-joon
% for Fitness data processing
% in my computer(com_name : scsi-PC)

clc; close all; clear all;

% Folder setting
Folder_path = 'E:\�л������� �м�\2017\rawdata\6. ��Ʈ�Ͻ�_���Ա��(180524 ����)\2017�� ���峻��\2017�� ���峻��';
cd(Folder_path)

% .xlsx reading
data = xlsread('2017_1��.xlsx');
% LUT = xlsread('LUT.xlsx');

% Cleaning(������ �����)
data(:,[2,3]) = [];
del = find(isnan(data(:,1)));
data(del,:) = [];

% Excel vlookup
% Origin_Stu_id = data(:,1);
% Stu_id = vlookup(Origin_Stu_id, LUT, 2); % �й��� ��ȯ��

Stu_id(Stu_id == 0) = NaN;
data(:,1) = Stu_id; % ��ȯ �ȵ� ���� NaN���� �ٲٰ� 1���� �߰�

del2 = find(isnan(data(:,1)));
data(del2,:) = [];
data(~isfinite(data))=0; % ��ȯ �ȵ� �� ��ü�� �����ϰ� 0���� �ٲ���

% Integrated data generation
k=0;
int = [];
for j=2:(size(data,2)-2)
    temp = data(:,[1 j]);
    for i=1:size(data,1)
        if temp(i,2) ~= 0
            k=k+1;
            int(k,1) = temp(i,1);
            int(k,2) = j-1;
            int(k,3) = temp(i,2);
        end
    end
end

for i=1:size(int,1)
    int(i, [4,5,6]) = [604,7,12];
end

xlswrite('FITNESS_12.xlsx',int);

function b=vlookup(m,lut,n)
[u,v]=size(m);
[w,~]=size(lut);
b=zeros(u,v);
    for i=1:u
        for j=1:v
            q=m(i,j);
            for k=1:w
                if q==lut(k)
                    b(i,j)=lut(k,n);
                    break
                end
            end
        end
    end
end