function [Ldistance] = p_dis(distance,num_can)
Ldistance = distance; % num_can * num_can ����

% ����ֱ���ҵ�����
for i = 1:num_can
    distance(i,i) = inf;
end

for i = 1:(num_can-1)
    for j = (i+1):num_can
        t = 1; % ��¼Ѱ�Ҵ���
        k(t) = j; % ��¼�ô��ҵĵ㣬���ӵ�j��ʼ��i��
        
%         id = num_can*(i-1) + i*(i+1)/2 + (j-1); % �ھ���������Ҷ�Ӧ�������
        turn = 1;
        dik = distance(i,k(t)); 
        dkk = distance(i,k(t));
        max_dkk = 0;
        
        while turn == 1
            
            % ���е���i�ľ��룬�ҳ�С�ڵ�i�� k(t)֮�����ĵ�
            line_t1 = find(distance(i,:) < dik); % ���ﲻ�����õ㱾��
            
            % �����������ĵ㼯Ϊ�գ��Ƕ���ij֮���L�������i����t����k(t)����
            if isempty(line_t1)
                if max_dkk == 0
                   Ldistance(i,j) = distance(i,j);
                else
                   Ldistance(i,j) = max_dkk;
                   Ldistance(j,i) = max_dkk;
                end
                break
            end
            % �ҳ�����������k(t)������С��һ��
            [~, id_temp] = min(distance(line_t1, k(t)));
            id_t = line_t1(id_temp);
            % ���������С�ĵ�Ϊk(t+1)
            t = t+1;
            k(t) = id_t;  
            dkk = distance(k(t),k(t-1));
            dik = distance(i,k(t));
            % ��ÿ�ε�������С��������ֵ�����¼����
            if dkk > max_dkk
%                 k(t), k(t-1);
                max_dkk = dkk;
%                 pause
            end
            
        end
        
    end
end


                
