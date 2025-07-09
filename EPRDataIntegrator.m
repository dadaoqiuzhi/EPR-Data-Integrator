%EPR/ESR谱图平滑和自由基浓度分析
%version 1, 2021.4.9
%刘强，核物理与化学研究所210室，forliubinqiang@163.com,15828636974

disp('每个条件的数据一起放进来，实现对一组的作差基线校正和数据平滑,基线校正是基于原始曲线所有数据')
close all;
n=input('请输入过滤级别n,推荐30合适:\n');
calmethod=input('\n请选择计算方法：1.直接一次和二次积分  2.峰高差乘以峰中心位置差:\n');
region=input('请输入积分或者计算区间,先小后大，多个区间连续输入即可，空格隔开,带引号：\n');
region=strtrim(region);region=strsplit(region);
if rem(length(region),2)~=0%mod()函数也可
    error('积分或者计算区间输入非法，并不是偶数，请检查！！！');
end
regiondata={};
for i=1:length(region)/2
    regiondata{i,1}=str2num(region{1,2*i-1});
    regiondata{i,2}=str2num(region{1,2*i});
end

%读入数据
data=xlsread('input_data.xlsx');
[row,col]=size(data);
%原始数据减去基线
for i=1:col/2
    data(:,2*i)=data(:,2*i)-mean(data(:,2*i));
end
%多图作图设置
datapro=[];
jj=1;figure(jj);k=1;
for i=1:col/2%平滑，矫正基线
    if ceil(i/6)>jj%控制为6小图
        jj=jj+1;
        figure(jj);
        k=1;
    end
    datapro(:,i)=medfilt1(data(:,2*i),n);
    datapro(:,i)=datapro(:,i)-mean(datapro(:,i));
    subplot(2,3,k)
    k=k+1;
    plot(data(:,1),data(:,2*i),'-k','linewidth',1.5);%原始数据
    hold on
    plot(data(:,1),datapro(:,i),'-r','linewidth',3);%平滑数据
    hold off
    legend( 'origin data','fitting data', 'Location', 'NorthEast' );
    xlabel( 'Field (G)' );ylabel( 'Intensity' )
    title('平滑和基线校正')
end


    

dataprocopy=datapro;curvesegment={};
for i=1:size(regiondata,1)%循环处理多峰谱图的每个区间
    datapro=dataprocopy;
    
    for j=1:row%截留积分区间
        if data(j,1)>=regiondata{i,1}
            regiondata{i,3}=j;
            break
        end
    end
    for k=1:row%截留积分区间
        if data(k,1)>=regiondata{i,2}
            regiondata{i,4}=row-k;
            break
        end
    end
    datafield=data(:,1);
    datapro(1:regiondata{i,3},:)=[];
    datafield(1:regiondata{i,3},:)=[];
    datapro(row-regiondata{i,3}-regiondata{i,4}+1:end,:)=[];
    datafield(row-regiondata{i,3}-regiondata{i,4}+1:end,:)=[];
    
    curvesegment{i,1}=datafield(:,1);
    for j=1:col/2 %存储每个区间的谱图好数据
        curvesegment{i,j+1}=datapro(:,j);
    end
end
%处理浓度
if calmethod==1
    integraI={};integraII={};
    for i=1:size(curvesegment,1)
        integraI{i,1}=curvesegment{i,1};
        integraII{i,1}=curvesegment{i,1};
        for j=1:1:size(curvesegment,2)-1
            integraI{i,j+1}=cumtrapz(integraI{i,1},curvesegment{i,j+1});
            integraII{i,j+1}=cumtrapz(integraII{i,1},integraI{i,j+1});
        end
    end
    
    for i=1:size(integraI,1)
        jj=jj+1;figure(jj);k=1;m=jj;%一次积分后画图，以同一区间的积分为一幅图
        for j=1:size(integraI,2)-1
            if ceil(j/6)>jj-m+1%控制为6小图
                jj=jj+1;
                figure(jj);
                k=1;
            end
            subplot(2,3,k)
            k=k+1;
            plot(integraI{i,1},integraI{i,j+1},'-k','linewidth',3);
            xlabel( 'Field (G)' );ylabel( 'Intensity-intI' )
            titlename=strcat('一次积分of数据',num2str(j),'(',num2str(i),',',num2str(j),')');
            title(titlename)
        end
    end
    
    integraIIdata_max=[];integraIIdata_end=[];
    for i=1:size(integraII,1)
        jj=jj+1;figure(jj);k=1;m=jj;%二次积分后画图，以同一区间的积分为一幅图
        for j=1:size(integraII,2)-1
            if ceil(j/6)>jj-m+1%控制为6小图
                jj=jj+1;
                figure(jj);
                k=1;
            end
            subplot(2,3,k)
            k=k+1;
            plot(integraII{i,1},integraII{i,j+1},'-k','linewidth',3);
            xlabel( 'Field (G)' );ylabel( 'Intensity-intII' )
            titlename=strcat('二次积分of数据',num2str(j),'(',num2str(i),',',num2str(j),')');
            title(titlename)
            fprintf('第%d组数据积分矩阵%s最大值为%f,积分最后值为%f\n',j,strcat('(',num2str(i),',',num2str(j),')'),max(integraII{i,j+1}),integraII{i,j+1}(end,1))
            integraIIdata_max(i,j)=max(integraII{i,j+1});
            integraIIdata_end(i,j)=integraII{i,j+1}(end,1);
        end
    end
    fprintf('\n一次积分和二次积分数据分别在integraI和integraII中，二次积分最大值和最后值储分别存在integraIIdata_max和integraIIdata_end中\n');
    clear calmethod col i j jj k m n row region regiondata
    
elseif calmethod==2
    concentration={};%存储峰高差和峰键距以及他们的乘积值
    concentrationdata=[];%只存储面积
    for i=1:size(curvesegment,1)
        for j=1:size(curvesegment,2)-1
            [height1,index1]=max(curvesegment{i,j+1}(:,1));
            [height2,index2]=min(curvesegment{i,j+1}(:,1));
            peak_height=height1-height2;
            peak_span=abs(curvesegment{i,1}(index2,1)-curvesegment{i,1}(index1,1));
            area=peak_height*peak_span;
            data_con=[peak_height,peak_span,area];
            concentration{i,j}=data_con(1,:);
            concentrationdata(i,j)=area;
        end
    end
    jj=jj+1;figure(jj);
    
    linesymb={'+','o','*','x','s','^','v','>','<','p','h'};
    linecolor={'r','g','b','c','m','y','k'};
    linetype={'-','--',':','-.'};
    
    for k=1:size(concentrationdata,1);
        rand1=round(10*rand())+1;rand2=round(6*rand())+1;rand3=round(3*rand())+1;%随机选择线型、颜色和符号
        linespecifer=strcat(linesymb{rand1},linecolor{rand2},linetype{rand3});
        plot(concentrationdata(k,:),linespecifer,'LineWidth',2);
        hold on
    end
    hold off
    xlabel( 'Sample No.' );ylabel( 'Intensity (a.u.)' )
    title('Concentration of radical/spin');
    legend_opt=size(regiondata,1);
    switch legend_opt
        case 1
            legend('区间1面积data');
        case 2
            legend('区间1面积data','区间2面积data');
        case 3
            legend('区间1面积data','区间2面积data','区间3面积data');
        case 4
            legend('区间1面积data','区间2面积data','区间3面积data','区间4面积data');
        case 5
            legend('区间1面积data','区间2面积data','区间3面积data','区间4面积data','区间5面积data');
        otherwise
            fprintf('\n计算区间超过5，不再显示标题，每条线是不同样品同一区间的计算面积！！！\n');
    end
    
    
    fprintf('基线校正和平滑处理后的数据储存在datapro中；分段数据储存在curvesegment中，峰高、峰间距和面积存储在concentration中，峰面积单独存储在concentrationdata中\n')
    clear area calmethod col data_con height1 height2 i j k m n index1 index2 peak_height peak_span row
end







