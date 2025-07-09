%EPR/ESR��ͼƽ�������ɻ�Ũ�ȷ���
%version 1, 2021.4.9
%��ǿ���������뻯ѧ�о���210�ң�forliubinqiang@163.com,15828636974

disp('ÿ������������һ��Ž�����ʵ�ֶ�һ����������У��������ƽ��,����У���ǻ���ԭʼ������������')
close all;
n=input('��������˼���n,�Ƽ�30����:\n');
calmethod=input('\n��ѡ����㷽����1.ֱ��һ�κͶ��λ���  2.��߲���Է�����λ�ò�:\n');
region=input('��������ֻ��߼�������,��С��󣬶�������������뼴�ɣ��ո����,�����ţ�\n');
region=strtrim(region);region=strsplit(region);
if rem(length(region),2)~=0%mod()����Ҳ��
    error('���ֻ��߼�����������Ƿ���������ż�������飡����');
end
regiondata={};
for i=1:length(region)/2
    regiondata{i,1}=str2num(region{1,2*i-1});
    regiondata{i,2}=str2num(region{1,2*i});
end

%��������
data=xlsread('input_data.xlsx');
[row,col]=size(data);
%ԭʼ���ݼ�ȥ����
for i=1:col/2
    data(:,2*i)=data(:,2*i)-mean(data(:,2*i));
end
%��ͼ��ͼ����
datapro=[];
jj=1;figure(jj);k=1;
for i=1:col/2%ƽ������������
    if ceil(i/6)>jj%����Ϊ6Сͼ
        jj=jj+1;
        figure(jj);
        k=1;
    end
    datapro(:,i)=medfilt1(data(:,2*i),n);
    datapro(:,i)=datapro(:,i)-mean(datapro(:,i));
    subplot(2,3,k)
    k=k+1;
    plot(data(:,1),data(:,2*i),'-k','linewidth',1.5);%ԭʼ����
    hold on
    plot(data(:,1),datapro(:,i),'-r','linewidth',3);%ƽ������
    hold off
    legend( 'origin data','fitting data', 'Location', 'NorthEast' );
    xlabel( 'Field (G)' );ylabel( 'Intensity' )
    title('ƽ���ͻ���У��')
end


    

dataprocopy=datapro;curvesegment={};
for i=1:size(regiondata,1)%ѭ����������ͼ��ÿ������
    datapro=dataprocopy;
    
    for j=1:row%������������
        if data(j,1)>=regiondata{i,1}
            regiondata{i,3}=j;
            break
        end
    end
    for k=1:row%������������
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
    for j=1:col/2 %�洢ÿ���������ͼ������
        curvesegment{i,j+1}=datapro(:,j);
    end
end
%����Ũ��
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
        jj=jj+1;figure(jj);k=1;m=jj;%һ�λ��ֺ�ͼ����ͬһ����Ļ���Ϊһ��ͼ
        for j=1:size(integraI,2)-1
            if ceil(j/6)>jj-m+1%����Ϊ6Сͼ
                jj=jj+1;
                figure(jj);
                k=1;
            end
            subplot(2,3,k)
            k=k+1;
            plot(integraI{i,1},integraI{i,j+1},'-k','linewidth',3);
            xlabel( 'Field (G)' );ylabel( 'Intensity-intI' )
            titlename=strcat('һ�λ���of����',num2str(j),'(',num2str(i),',',num2str(j),')');
            title(titlename)
        end
    end
    
    integraIIdata_max=[];integraIIdata_end=[];
    for i=1:size(integraII,1)
        jj=jj+1;figure(jj);k=1;m=jj;%���λ��ֺ�ͼ����ͬһ����Ļ���Ϊһ��ͼ
        for j=1:size(integraII,2)-1
            if ceil(j/6)>jj-m+1%����Ϊ6Сͼ
                jj=jj+1;
                figure(jj);
                k=1;
            end
            subplot(2,3,k)
            k=k+1;
            plot(integraII{i,1},integraII{i,j+1},'-k','linewidth',3);
            xlabel( 'Field (G)' );ylabel( 'Intensity-intII' )
            titlename=strcat('���λ���of����',num2str(j),'(',num2str(i),',',num2str(j),')');
            title(titlename)
            fprintf('��%d�����ݻ��־���%s���ֵΪ%f,�������ֵΪ%f\n',j,strcat('(',num2str(i),',',num2str(j),')'),max(integraII{i,j+1}),integraII{i,j+1}(end,1))
            integraIIdata_max(i,j)=max(integraII{i,j+1});
            integraIIdata_end(i,j)=integraII{i,j+1}(end,1);
        end
    end
    fprintf('\nһ�λ��ֺͶ��λ������ݷֱ���integraI��integraII�У����λ������ֵ�����ֵ���ֱ����integraIIdata_max��integraIIdata_end��\n');
    clear calmethod col i j jj k m n row region regiondata
    
elseif calmethod==2
    concentration={};%�洢��߲�ͷ�����Լ����ǵĳ˻�ֵ
    concentrationdata=[];%ֻ�洢���
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
        rand1=round(10*rand())+1;rand2=round(6*rand())+1;rand3=round(3*rand())+1;%���ѡ�����͡���ɫ�ͷ���
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
            legend('����1���data');
        case 2
            legend('����1���data','����2���data');
        case 3
            legend('����1���data','����2���data','����3���data');
        case 4
            legend('����1���data','����2���data','����3���data','����4���data');
        case 5
            legend('����1���data','����2���data','����3���data','����4���data','����5���data');
        otherwise
            fprintf('\n�������䳬��5��������ʾ���⣬ÿ�����ǲ�ͬ��Ʒͬһ����ļ������������\n');
    end
    
    
    fprintf('����У����ƽ�����������ݴ�����datapro�У��ֶ����ݴ�����curvesegment�У���ߡ����������洢��concentration�У�����������洢��concentrationdata��\n')
    clear area calmethod col data_con height1 height2 i j k m n index1 index2 peak_height peak_span row
end







