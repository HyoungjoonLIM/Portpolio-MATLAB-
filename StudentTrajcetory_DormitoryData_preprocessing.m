
%% For unstructured data + structured dta
% Folder setting
Folder_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사\출입기록(new)\C동';
cd(Folder_path)
Folder_ls = ls;
Folder_ls(1:2,:) = [];

% Folder exploring
Result=[];
cnt = 0;
prb = 1;
for i = 1:size(Folder_ls,1)
    i
    cd([Folder_path,'\',Folder_ls(i,:)])
    File_ls1 = ls;
    File_ls1(1:2,:) = [];
    
    % File exploring
    for j = 1:size(File_ls1,1)
            j
            [num, txt, raw] = xlsread(File_ls1(j,:));
            % Unstructred data 1
                raw(1:5,:) = [];

                ind_1 = find(strcmp(raw(1,:),'일자'));
                ind_2 = find(strcmp(raw(1,:),'시간'));
                ind_3 = find(strcmp(raw(1,:),'건물'));
                ind_4 = find(strcmp(raw(1,:),'출입문'));
                ind_5 = find(strcmp(raw(1,:),'사원명'));
                ind_6 = find(strcmp(raw(1,:),'카드번호'));

                raw = cat(2,raw(:,ind_1),raw(:,ind_2),raw(:,ind_3),raw(:,ind_4),raw(:,ind_5),raw(:,ind_6));

                raw(1,:) = [];
                              
                % Find rows: '사원명' = NaN
                nanrow = [];
                for k = 1:size(raw,1)
                    if isnan(raw{k,5}) == 1
                        nanrow = [nanrow;k];
                    end
                end
                
                raw(nanrow,:)=[];
                
                % Excel exploring
                for n = 1:size(raw,1)
                    if ischar(raw{n,1}) == 0
                        raw{n-1,5} = strcat(raw{n-1,5},raw{n,5});
                    end
                end
                
                % Find rows: '사원명' = NaN
                nanrow2 = [];
                for m = 1:size(raw,1)
                    if isnan(raw{m,1}) == 1
                        nanrow2 = [nanrow2;m];
                    end
                end
                raw(nanrow2,:)=[];
                
            Result = [Result;raw];
            clear raw
            clear num
            clear txt
        
    end
    clear File_ls
end

for i = 1:size(Result,1)
    i
    logic(i) = length(Result{i,5}) > 1;
end

logic = logic';

Result_C = cell2table(Result);

clear Result

% save_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사';
% cd(save_path)
% csvwrite('dorm_B.csv', Result_B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For unstructured data + structured dta
% Folder setting
Folder_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사\출입기록(new)\D동 주출입구';
cd(Folder_path)
Folder_ls = ls;
Folder_ls(1:2,:) = [];

% Folder exploring
Result=[];
cnt = 0;
prb = 1;
for i = 1:size(Folder_ls,1)
    i
    cd([Folder_path,'\',Folder_ls(i,:)])
    File_ls1 = ls;
    File_ls1(1:2,:) = [];
    
    % File exploring
    for j = 1:size(File_ls1,1)
            j
            [num, txt, raw] = xlsread(File_ls1(j,:));
            % Unstructred data 1
                raw(1:5,:) = [];

                ind_1 = find(strcmp(raw(1,:),'일자'));
                ind_2 = find(strcmp(raw(1,:),'시간'));
                ind_3 = find(strcmp(raw(1,:),'건물'));
                ind_4 = find(strcmp(raw(1,:),'출입문'));
                ind_5 = find(strcmp(raw(1,:),'사원명'));
                ind_6 = find(strcmp(raw(1,:),'카드번호'));

                raw = cat(2,raw(:,ind_1),raw(:,ind_2),raw(:,ind_3),raw(:,ind_4),raw(:,ind_5),raw(:,ind_6));

                raw(1,:) = [];
                              
                % Find rows: '사원명' = NaN
                nanrow = [];
                for k = 1:size(raw,1)
                    if isnan(raw{k,5}) == 1
                        nanrow = [nanrow;k];
                    end
                end
                
                raw(nanrow,:)=[];
                
                % Excel exploring
                for n = 1:size(raw,1)
                    if ischar(raw{n,1}) == 0
                        raw{n-1,5} = strcat(raw{n-1,5},raw{n,5});
                    end
                end
                
                % Find rows: '사원명' = NaN
                nanrow2 = [];
                for m = 1:size(raw,1)
                    if isnan(raw{m,1}) == 1
                        nanrow2 = [nanrow2;m];
                    end
                end
                raw(nanrow2,:)=[];
                
            Result = [Result;raw];
            clear raw
            clear num
            clear txt
        
    end
    clear File_ls
end

for i = 1:size(Result,1)
    i
    logic(i) = length(Result{i,5}) > 1;
end

logic = logic';

Result_D_main = cell2table(Result);

clear Result

% save_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사';
% cd(save_path)
% csvwrite('dorm_B.csv', Result_B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For unstructured data + structured dta
% Folder setting
Folder_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사\출입기록(new)\D동 부출입구';
cd(Folder_path)
Folder_ls = ls;
Folder_ls(1:2,:) = [];

% Folder exploring
Result=[];
cnt = 0;
prb = 1;
for i = 1:size(Folder_ls,1)
    i
    cd([Folder_path,'\',Folder_ls(i,:)])
    File_ls1 = ls;
    File_ls1(1:2,:) = [];
    
    % File exploring
    for j = 1:size(File_ls1,1)
            j
            [num, txt, raw] = xlsread(File_ls1(j,:));
            % Unstructred data 1
                raw(1:5,:) = [];

                ind_1 = find(strcmp(raw(1,:),'일자'));
                ind_2 = find(strcmp(raw(1,:),'시간'));
                ind_3 = find(strcmp(raw(1,:),'건물'));
                ind_4 = find(strcmp(raw(1,:),'출입문'));
                ind_5 = find(strcmp(raw(1,:),'사원명'));
                ind_6 = find(strcmp(raw(1,:),'카드번호'));

                raw = cat(2,raw(:,ind_1),raw(:,ind_2),raw(:,ind_3),raw(:,ind_4),raw(:,ind_5),raw(:,ind_6));

                raw(1,:) = [];
                              
                % Find rows: '사원명' = NaN
                nanrow = [];
                for k = 1:size(raw,1)
                    if isnan(raw{k,5}) == 1
                        nanrow = [nanrow;k];
                    end
                end
                
                raw(nanrow,:)=[];
                
                % Excel exploring
                for n = 1:size(raw,1)
                    if ischar(raw{n,1}) == 0
                        raw{n-1,5} = strcat(raw{n-1,5},raw{n,5});
                    end
                end
                
                % Find rows: '사원명' = NaN
                nanrow2 = [];
                for m = 1:size(raw,1)
                    if isnan(raw{m,1}) == 1
                        nanrow2 = [nanrow2;m];
                    end
                end
                raw(nanrow2,:)=[];
                
            Result = [Result;raw];
            clear raw
            clear num
            clear txt
        
    end
    clear File_ls
end

for i = 1:size(Result,1)
    i
    logic(i) = length(Result{i,5}) > 1;
end

logic = logic';

Result_D_sub = cell2table(Result);

clear Result

% save_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사';
% cd(save_path)
% csvwrite('dorm_B.csv', Result_B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For unstructured data + structured dta
% Folder setting
Folder_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사\출입기록(new)\E동';
cd(Folder_path)
Folder_ls = ls;
Folder_ls(1:2,:) = [];

% Folder exploring
Result=[];
cnt = 0;
prb = 1;
for i = 1:size(Folder_ls,1)
    i
    cd([Folder_path,'\',Folder_ls(i,:)])
    File_ls1 = ls;
    File_ls1(1:2,:) = [];
    
    % File exploring
    for j = 1:size(File_ls1,1)
            j
            [num, txt, raw] = xlsread(File_ls1(j,:));
            % Unstructred data 1
                raw(1:5,:) = [];

                ind_1 = find(strcmp(raw(1,:),'일자'));
                ind_2 = find(strcmp(raw(1,:),'시간'));
                ind_3 = find(strcmp(raw(1,:),'건물'));
                ind_4 = find(strcmp(raw(1,:),'출입문'));
                ind_5 = find(strcmp(raw(1,:),'사원명'));
                ind_6 = find(strcmp(raw(1,:),'카드번호'));

                raw = cat(2,raw(:,ind_1),raw(:,ind_2),raw(:,ind_3),raw(:,ind_4),raw(:,ind_5),raw(:,ind_6));

                raw(1,:) = [];
                              
                % Find rows: '사원명' = NaN
                nanrow = [];
                for k = 1:size(raw,1)
                    if isnan(raw{k,5}) == 1
                        nanrow = [nanrow;k];
                    end
                end
                
                raw(nanrow,:)=[];
                
                % Excel exploring
                for n = 1:size(raw,1)
                    if ischar(raw{n,1}) == 0
                        raw{n-1,5} = strcat(raw{n-1,5},raw{n,5});
                    end
                end
                
                % Find rows: '사원명' = NaN
                nanrow2 = [];
                for m = 1:size(raw,1)
                    if isnan(raw{m,1}) == 1
                        nanrow2 = [nanrow2;m];
                    end
                end
                raw(nanrow2,:)=[];
                
            Result = [Result;raw];
            clear raw
            clear num
            clear txt
        
    end
    clear File_ls
end

for i = 1:size(Result,1)
    i
    logic(i) = length(Result{i,5}) > 1;
end

logic = logic';

Result_E = cell2table(Result);

clear Result

% save_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사';
% cd(save_path)
% csvwrite('dorm_B.csv', Result_B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For unstructured data + structured dta
% Folder setting
Folder_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사\출입기록(new)\F동';
cd(Folder_path)
Folder_ls = ls;
Folder_ls(1:2,:) = [];

% Folder exploring
Result=[];
cnt = 0;
prb = 1;
for i = 1:size(Folder_ls,1)
    i
    cd([Folder_path,'\',Folder_ls(i,:)])
    File_ls1 = ls;
    File_ls1(1:2,:) = [];
    
    % File exploring
    for j = 1:size(File_ls1,1)
            j
            [num, txt, raw] = xlsread(File_ls1(j,:));
            % Unstructred data 1
                raw(1:5,:) = [];

                ind_1 = find(strcmp(raw(1,:),'일자'));
                ind_2 = find(strcmp(raw(1,:),'시간'));
                ind_3 = find(strcmp(raw(1,:),'건물'));
                ind_4 = find(strcmp(raw(1,:),'출입문'));
                ind_5 = find(strcmp(raw(1,:),'사원명'));
                ind_6 = find(strcmp(raw(1,:),'카드번호'));

                raw = cat(2,raw(:,ind_1),raw(:,ind_2),raw(:,ind_3),raw(:,ind_4),raw(:,ind_5),raw(:,ind_6));

                raw(1,:) = [];
                              
                % Find rows: '사원명' = NaN
                nanrow = [];
                for k = 1:size(raw,1)
                    if isnan(raw{k,5}) == 1
                        nanrow = [nanrow;k];
                    end
                end
                
                raw(nanrow,:)=[];
                
                % Excel exploring
                for n = 1:size(raw,1)
                    if ischar(raw{n,1}) == 0
                        raw{n-1,5} = strcat(raw{n-1,5},raw{n,5});
                    end
                end
                
                % Find rows: '사원명' = NaN
                nanrow2 = [];
                for m = 1:size(raw,1)
                    if isnan(raw{m,1}) == 1
                        nanrow2 = [nanrow2;m];
                    end
                end
                raw(nanrow2,:)=[];
                
            Result = [Result;raw];
            clear raw
            clear num
            clear txt
        
    end
    clear File_ls
end

for i = 1:size(Result,1)
    i
    logic(i) = length(Result{i,5}) > 1;
end

logic = logic';

Result_F = cell2table(Result);

clear Result

% save_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사';
% cd(save_path)
% csvwrite('dorm_B.csv', Result_B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For unstructured data + structured dta
% Folder setting
Folder_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사\출입기록(new)\G동';
cd(Folder_path)
Folder_ls = ls;
Folder_ls(1:2,:) = [];

% Folder exploring
Result=[];
cnt = 0;
prb = 1;
for i = 1:size(Folder_ls,1)
    i
    cd([Folder_path,'\',Folder_ls(i,:)])
    File_ls1 = ls;
    File_ls1(1:2,:) = [];
    
    % File exploring
    for j = 1:size(File_ls1,1)
            j
            [num, txt, raw] = xlsread(File_ls1(j,:));
            % Unstructred data 1
                raw(1:5,:) = [];

                ind_1 = find(strcmp(raw(1,:),'일자'));
                ind_2 = find(strcmp(raw(1,:),'시간'));
                ind_3 = find(strcmp(raw(1,:),'건물'));
                ind_4 = find(strcmp(raw(1,:),'출입문'));
                ind_5 = find(strcmp(raw(1,:),'사원명'));
                ind_6 = find(strcmp(raw(1,:),'카드번호'));

                raw = cat(2,raw(:,ind_1),raw(:,ind_2),raw(:,ind_3),raw(:,ind_4),raw(:,ind_5),raw(:,ind_6));

                raw(1,:) = [];
                              
                % Find rows: '사원명' = NaN
                nanrow = [];
                for k = 1:size(raw,1)
                    if isnan(raw{k,5}) == 1
                        nanrow = [nanrow;k];
                    end
                end
                
                raw(nanrow,:)=[];
                
                % Excel exploring
                for n = 1:size(raw,1)
                    if ischar(raw{n,1}) == 0
                        raw{n-1,5} = strcat(raw{n-1,5},raw{n,5});
                    end
                end
                
                % Find rows: '사원명' = NaN
                nanrow2 = [];
                for m = 1:size(raw,1)
                    if isnan(raw{m,1}) == 1
                        nanrow2 = [nanrow2;m];
                    end
                end
                raw(nanrow2,:)=[];
                
            Result = [Result;raw];
            clear raw
            clear num
            clear txt
        
    end
    clear File_ls
end

for i = 1:size(Result,1)
    i
    logic(i) = length(Result{i,5}) > 1;
end

logic = logic';

Result_G = cell2table(Result);

clear Result


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% csvwrite
save_path = 'E:\학생데이터 분석\2017\rawdata\5. 기숙사\기숙사';
cd(save_path)
writetable(Result_A, 'dorm_A.csv');
writetable(Result_B, 'dorm_B.csv');
writetable(Result_C, 'dorm_C.csv');
writetable(Result_D_main, 'dorm_D1.csv');
writetable(Result_D_sub, 'dorm_D2.csv');
writetable(Result_E, 'dorm_E.csv');
writetable(Result_F, 'dorm_F.csv');
writetable(Result_G, 'dorm_G.csv');





