clear all;  
clc;
close all;
load A2.txt;

tic;
data = A2(:,1:(end-1));

%% centralize and scale the data
data = data - repmat(mean(data),size(data,1),1);
data = data/max(max(abs(data)));

true_label = A2(:,end);

halo = true_label;
NClusters = max(halo);
% twoA2 = [[two_A2(xx1,[1,2]),yy1];[two_A2(xx2,[1,2]),-yy2]];

figure
hold on
cmap=colormap;
for i=1:NClusters
   ic=int8((i*64.)/(NClusters*1.));

   if NClusters<=12
       switch i
            case 1
                  plot(data(halo==i,1),data(halo==i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 2
                  plot(data(halo==i,1),data(halo==i,2),'*','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 3
                  plot(data(halo==i,1),data(halo==i,2),'x','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 4
                  plot(data(halo==i,1),data(halo==i,2),'+','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 5
                  plot(data(halo==i,1),data(halo==i,2),'s','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 6
                  plot(data(halo==i,1),data(halo==i,2),'d','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 7
                  plot(data(halo==i,1),data(halo==i,2),'v','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 8
                  plot(data(halo==i,1),data(halo==i,2),'^','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 9
                  plot(data(halo==i,1),data(halo==i,2),'<','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 10
                  plot(data(halo==i,1),data(halo==i,2),'>','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 11
                  plot(data(halo==i,1),data(halo==i,2),'p','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 12
                  plot(data(halo==i,1),data(halo==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
       end

    else
      plot(data(halo==i,1),data(halo==i,2),'.','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));	
    end % end of if 					  

end %end of for
set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.28]);
set(gca,'FontName','Times New Roman','FontSize',14)
title('Local Center','FontSize',16)	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = size(data,2);
num = size(data,1);

% build kd-tree
NS = createns(data,'NSMethod','kdtree');

% search knn for every points
knn = 20 + floor((num/1000)*5)
find_k = 19;
% max_time = 2; % 用以判断是否为邻近群
% max_nei = 7;
num_k = knn + 1 ;  %knn+1,for thr first con
nei = zeros(num, num_k); % record nei index
distance = zeros(num, num_k); % record nei distanceance
% min_ther = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 初次筛选候选中心点
% 遍历，对每个点求k近邻以及k近邻与其的平均距离，代表密度
for i = 1:num
    [nei(i,:), distance(i,:)] = knnsearch(NS,data(i,:),'k',num_k);
    mean_distance(i,:) = mean(distance(i,:));  %% 这里的平均距离的倒数，定义为rho
end
nei(:,1) = []; % 消除该点本身
distance(:,1) = [];

% 对密度排序，从密度大的向下遍历，确定局部中心
[mean_distance_sort, id_all_sort] = sort(mean_distance, 'ascend'); % 所有点根据平均距离排序

center_id1 = [];
for i = 1: num 
    p1 = id_all_sort(i);
    distance_p1 = mean_distance_sort(i);
    same_point = intersect(nei(p1,1:find_k),center_id1);
    distance_nei = mean_distance(nei(p1,1:find_k));
   
     % 如果该点比起大部分邻居点的密度都小，则不作为can
    if isempty (same_point) && (distance_p1 <= mean(distance_nei))
        center_id1 = [center_id1 ,p1];
    end

end   

center_can = center_id1';
num_can = length(center_can);

% 画骨架点
if dim == 2
    id = center_can;
    for i = 1:length(id)
        i;
        plot(data(id(i),1),data(id(i),2), '*','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','k');
        plot(data(id(i),1),data(id(i),2), 'o','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');
    %     pause
    end
end

%%%%%%%%%%%%%%% CK %%%%%%%%%%%%%%%%%%%%%%%%%
% center_can = id_all_sort; % change to CK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_can = data(center_can,:);
num_can = size(data_can,1);

%% 用can计算density peaks               
mdistance = [];
% 欧式距离矩阵				
for i=1:num_can  
      mdistance_i = [];
      for j=i+1:num_can 
            mdistance_i = [mdistance_i;[i,j, sqrt(sum((data_can(i,:) - data_can(j,:)) .^ 2)) ]];
      end
      mdistance = [mdistance; mdistance_i];
end


%% compute and save diatance matrix if dont exists			
save A2DPCmdistance.distance   mdistance -ascii                
                
%% load diatance matrix if exists				
load A2DPCmdistance.distance;	
xx = A2DPCmdistance;

ND=max(xx(:,2));
NL=max(xx(:,1));				  
				            
if (NL>ND)
    ND=NL;  %% num of points
end

N=size(xx,1); %% ND*(ND-1)/2,the number of pair distanceances
			 		
%% matrix of distanceance?ND*ND
distance = zeros(ND,ND);

%% matrix of distanceance?ND*ND, symmetry
for i=1:N
      ii=xx(i,1);
      jj=xx(i,2);
      distance(ii,jj)=xx(i,3);
      distance(jj,ii)=xx(i,3);
end

%% 计算新联通距离
[Ldistance] = p_dis(distance,num_can);

%% 设计新距离（结合Ldis 和 Edis）

% Cdistance = distance; % change to  DPC
Cdistance = Ldistance;

%% initialize?rho: density of points 
for i=1:ND
    mean_distance_can(i) = mean_distance(center_can(i));
    rho(i) = 1/mean_distance_can(i);
end

%% find the max distanceance?
maxd=max(max(Cdistance));

%% rank rho by descend?
[rho_sorted,ordrho] = sort(rho,'descend');
rho_or = rho_sorted;
ordrho_or = ordrho;

%% deal with point with max rho?
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0; 

%% compute the delta(relative distanceance), find nneigh for points
for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(Cdistance(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii)) = Cdistance(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
       % nneigh record the point with higher rho and most close to ordrho(ii)
     end
   end
end

%% give max rho point max delta 
delta(ordrho(1))=max(delta(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
%% decision graph 

disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')

fid = fopen('DECISION_GRAPH', 'w');
for i=1:ND
   fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
end

% Select a rectangle enclosing cluster centers
disp('Select a rectangle enclosing cluster centers')

scrsz = get(0,'ScreenSize');			 

fig = figure;
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.36]);
set(gca,'FontName','Times New Roman','FontSize',14)
title ('Decision Graph of A2','FontSize',17.0)
box off;
xlabel ('\rho','FontSize',15.0)
ylabel ('\delta','FontSize',15.0)

rect = getrect(fig);
rhomin=rect(1);
deltamin=rect(2); 

% initialize number of cluster
NCLUST=0;

%%%%%%%%%%%%%%%%%%%%%%% cluster center%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

%% cl(i)=j: point(i) belong to cluster j, initialize cl = -1;
for i=1:ND
  cl(i)=-1;
end

% find cluster center?
for i=1:ND
  if ( (rho(i)>rhomin) && (delta(i)>deltamin))
     NCLUST=NCLUST+1;
     cl(i)=NCLUST; %% cl(i) is center
     icl(NCLUST)=i;%% icl record index of cluster center
  end
end

fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
a = axis;

% limits of current figure
xmin=a(1);
xmax=a(2);
ymin=a(3);
ymax=a(4);

% makes grid
options.gridx = 50;
options.gridy = 50;
[X,Y] = meshgrid(xmin:(xmax-xmin)/options.gridx:xmax,...
                 ymin:(ymax-ymin)/options.gridy:ymax);

% make testing patterns covering whole grid
tst_data=[reshape(X',1,prod(size(X)));reshape(Y',1,prod(size(Y)))];
dec_fun= tst_data(1,:).*tst_data(2,:);
% reshape dec_fun
Z = reshape(dec_fun,size(X,1),size(X,2))';
% smooth shading
hold on
%contour(X,Y,Z,1,'k');

% draw cluster centers with different color
cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
   hold on
   cmap(ic,:)
   contour(X,Y,Z,[rho(icl(i))*delta(icl(i)) rho(icl(i))*delta(icl(i)) rho(icl(i))*delta(icl(i))],'linecolor',cmap(ic,:));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
disp('Performing assignation')

        %% clusting non_center points, traverse by rho?
for i=1:ND
      if (cl(ordrho(i))==-1)
         cl(ordrho(i))=cl(nneigh(ordrho(i)));
      end
end		

halo_can = cl';

%% for every cluster
for i=1:NCLUST
      nc=0; %% the num of all points in cluster
      nh=0; %% the num of core points in cluster 
      for j=1:ND
            if (cl(j)==i)
                 nc=nc+1;
            end
            if (halo_can(j)==i) % non_outlier
                 nh=nh+1;
            end
      end
      fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of DPC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ic=0;
if (dim==2)
    figure
    hold on
    cmap=colormap;
    for i=1:NCLUST
       ic=int8((i*64.)/(NCLUST*1.));
       if NCLUST<=12
           switch i
            case 1
                  plot(data_can(icl(i),1),data_can(icl(i),2),'h','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 2
                  plot(data_can(icl(i),1),data_can(icl(i),2),'h','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'*','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 3
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'x','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 4
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'+','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 5
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'s','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 6
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'d','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 7
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'v','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 8
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'^','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 9
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'<','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 10
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'>','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 11
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'p','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 12
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            otherwise
                  plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
           end

      else
           plot(data_can(icl(i),1),data_can(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
           plot(data_can(halo_can==i,1),data_can(halo_can==i,2),'.','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));	
      end 					  

    end
    plot(data_can(halo_can==0,1),data_can(halo_can==0,2),'.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');	

    id = icl;
    for i = 1:length(id)
    %     i
        plot(data_can(id(i),1),data_can(id(i),2), '*','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');
        plot(data_can(id(i),1),data_can(id(i),2), 'o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
    %     pause
    end
    set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.28]);
    set(gca,'FontName','Times New Roman','FontSize',14)
    title('A2','FontSize',16.0)

end



%% 分配剩余的点
% 将can的标签先加入总体
halo = zeros(num,1) - 1;

for i = 1 : num_can
    halo(center_can(i)) = halo_can(i);
    for j = 1:(find_k*2)
        if halo(nei(center_can(i),j)) == -1
            halo(nei(center_can(i),j))= halo_can(i);
        end
    end
end

% 将剩余点分配给距离最近的有标签点
for i = 1:num
    ii = id_all_sort(i);
    if  halo(ii) == -1
        for j = 1:knn
            if halo(nei(ii,j)) > 0
                halo(ii) = halo(nei(ii,j));
                break
            end
        end
    end
end


for i = 1:length(icl)
    icl(i) = center_can(icl(i));
end

if (dim==2)
    figure
    hold on
    cmap=colormap;
    for i=1:NCLUST
       ic=int8((i*64.)/(NCLUST*1.));
       if NCLUST<=12
           switch i
            case 1
                  plot(data(icl(i),1),data(icl(i),2),'h','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 2
                  plot(data(icl(i),1),data(icl(i),2),'h','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'*','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 3
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'x','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 4
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'+','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 5
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'s','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 6
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'d','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 7
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'v','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 8
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'^','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 9
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'<','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 10
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'>','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 11
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'p','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            case 12
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            otherwise
                  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                  plot(data(halo==i,1),data(halo==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
           end

      else
           plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
           plot(data(halo==i,1),data(halo==i,2),'.','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));	
      end 					  

  end
  plot(data(halo==0,1),data(halo==0,2),'.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');	

    hold on
    id = icl;
    for i = 1:length(id)
    %     i
        plot(data(id(i),1),data(id(i),2), '*','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');
        plot(data(id(i),1),data(id(i),2), 'o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
    %     pause
    end
    set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.28]);
    set(gca,'FontName','Times New Roman','FontSize',14)
    title('A2','FontSize',16.0) 
  
end

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FMI,ARI,NMI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    DPCFMI = FMI(true_label(:),halo(:));
    fprintf('FMI value of DPC on A2 dataset %i \n', DPCFMI);
    [cluster_acc,randindex,DPCARI] = ARI(true_label(:),halo(:));
    fprintf('ARI value of DPC on A2 dataset %i \n', DPCARI);
    DPCNMI = NMI(true_label(:),halo(:));
    fprintf('NMI value of DPC on A2 dataset %i \n', DPCNMI);
						 
				