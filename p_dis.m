function [Ldistance] = p_dis(distance,num_can)
Ldistance = distance; % num_can * num_can 矩阵

% 避免直接找到自身
for i = 1:num_can
    distance(i,i) = inf;
end

for i = 1:(num_can-1)
    for j = (i+1):num_can
        t = 1; % 记录寻找次数
        k(t) = j; % 记录该次找的点，即从点j开始向i找
        
%         id = num_can*(i-1) + i*(i+1)/2 + (j-1); % 在距离矩阵中找对应距离的行
        turn = 1;
        dik = distance(i,k(t)); 
        dkk = distance(i,k(t));
        max_dkk = 0;
        line_t2 = []; 
        while turn == 1
            
            % 所有点与i的距离，找出小于点i与 k(t)之间距离的点
            line_t1 = find(distance(i,~line_t2) < dik); % 这里不包括该点本身
            line_t2_temp = find(distance(i,:) > dik);
            line_t2 = [line_t2, line_t2_temp];
            % 若符合条件的点集为空，那定义ij之间的L距离就是i点与t代的k(t)距离
            if isempty(line_t1)
                if max_dkk == 0
                   Ldistance(i,j) = distance(i,j);
                else
                   Ldistance(i,j) = max_dkk;
                   Ldistance(j,i) = max_dkk;
                end
                break
            end
            % 找出上述点中与k(t)距离最小的一个
            [~, id_temp] = min(distance(line_t1, k(t)));
            id_t = line_t1(id_temp);
            % 定义距离最小的点为k(t+1)
            t = t+1;
            k(t) = id_t;  
            dkk = distance(k(t),k(t-1));
            dik = distance(i,k(t));
            % 将每次迭代中最小距离的最大值距离记录下来
            if dkk > max_dkk
%                 k(t), k(t-1);
                max_dkk = dkk;
%                 pause
            end
            
        end
        
    end
end


                
