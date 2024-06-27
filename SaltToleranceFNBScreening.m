
function SaltToleranceFNBScreening()
%--------------------------------------------------------------------------
% First Screening
traitPath='T1d-3d-5d-SpectRGBFeatures84.mat';
data=load(traitPath);
FBNTrait=data.AllTrait; % 1080*84*3 (FNB:1-1044  CK:1045-1062  SS:1063-1080)

nSamlpe=size(FBNTrait,1);
nA17=36;
AllResult=nan(nSamlpe-nA17,4);
for m=1:3 %
     [TBT_Class_Result,wt,id,Bb]=Class_TBT_Select(FBNTrait(:,:,m));
     AllResult(:,m+1)=TBT_Class_Result;
end
AllResult(AllResult==2)=0; %unselected FNB set 0
id=1:1044;
AllResult(:,1)=id';
SumClass=sum(AllResult(:,2:4),2);
FirstSelectedFNBID=AllResult(find(SumClass>2),1); % the id number of selected Salt-tolerant FBN 
clear traitPath data FBNTrait nA17 AllResult SumClass
%% ========================================================================
% Second Screening
metaPath='MetaTBT_5d-MSR_101Meta.mat';
data=load(metaPath);
FBNTmeta=data.Meta101; 

SeletID=zeros(nSamlpe,1);
SeletID(FirstSelectedFNBID,1)=1; 
SeletID(1045:1080,1)=1;  % A17CK  A17SS


ClassMateData=FBNTmeta(logical(SeletID),:);
[TBT_Class_Result,wt,id,Bb]=Class_TBT_Select(ClassMateData);
AllResult(:,2)=double(TBT_Class_Result);
AllResult(AllResult==2)=0;
AllResult(:,1)=FirstSelectedFNBID;
SecondSelectedFNBID=AllResult(find(AllResult(:,2)>0),1); 
clear metaPath data SeletID ClassMateData AllResult
%% =========================================================================
% PCA image of A17CK/A17SS/Salt-toloreant FNB/ unselected FNB
id=zeros(nSamlpe,1);
StFNBid=id;
StFNBid(SecondSelectedFNBID,1)=1; % logical value of salt-tolerant FNB id 
UstFNBid=~StFNBid;         
UstFNBid(1045:1080,1)=0;         % logical value of unselected FNB id
A17CK=id;
A17CK(1045:1062,1)=1;             % logical value of A17CK id

A17SS=id;
A17SS(1063:1080,1)=1;             % logical value of A17SS id


NorMate = normlization(FBNTmeta,1);  %标准化
[coeff,score,latent,tsquared,explained] = pca(NorMate);
%--------------------------------------------------------------------------
figure(8);
set(gcf,'Units','centimeter','Position',[5 0 18 20]);%2H 3L

alpha=0.05;
distribution='norm';
MarkerSize=6;
scatter(0,0,0.01,'filled','o','MarkerFaceColor',[0.30,0.75,0.93],'MarkerEdgeColor',[0 0 0]) % unselected FNB
hold on 
scatter(0,0,0.01,'filled','^','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0 0 0])       % salt-tolerant FNB
hold on
scatter(0,0,0.01,'filled','square','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0 0])     % A17CK
hold on 
scatter(0,0,0.01,'filled','diamond','MarkerFaceColor',[0.72,0.27,1.00],'MarkerEdgeColor',[0 0 1]); % A17SS

hold on
B3=plot([0,0],[-20,20],'--','LineWidth',1.1,'color',[0,0,0]);
set(B3,'handlevisibility','off');
hold on
B3=plot([-20,20],[0,0],'--','LineWidth',1.1,'color',[0,0,0]);
set(B3,'handlevisibility','off');


% hold on
hold on
P6_ConfidenceRegionFunction(score(logical(UstFNBid),1:2),alpha,distribution,[0.30,0.75,0.93],[0 0 0],'o',20);
hold on
P6_ConfidenceRegionFunction(score(logical(StFNBid),1:2),alpha,distribution,[1 0 0],[0 0 0],'^',25);
hold on
P6_ConfidenceRegionFunction(score(logical(A17CK),1:2),alpha,distribution,[0 1 0],[0 0 0],'square',25);
hold on
P6_ConfidenceRegionFunction(score(logical(A17SS),1:2),alpha,distribution,[0.72,0.27,1.00],[0 0 1],'diamond',35);

legend({'Unselected mutants','Salt-tolerant mutants','A17-CK','A17-SS'},'box','off',...
     'Position',[0.3 0.1 0.41 0.02],...
    'Orientation','horizontal');xlabel('Principal Component 1')
str1=sprintf('PC1(%4.2f%%)',explained(1,1));
str2=sprintf('PC2(%4.2f%%)',explained(2,1));
xlabel(str1)
ylabel(str2)
% xlabel('PC1(35.18%)')
% ylabel('PC2(19.88%)')
box on
axis equal %square %
ytick=[-15,0,15];
xtick=[-15,0,15];

set(gca,'yLim',[-20,15],'xTick',xtick,'XLim',[-20,15],'yTick',ytick,'Fontname', 'Arial',...
    'Fontsize', 10,'FontWeight','bold');
disp(0);

%% ========================================================================
% PCA3D display of Selected Salt-tolerant FNB
PCA_TBTdata=score(SecondSelectedFNBID,1:3);%1:3
%--------------------------------------------------------------------------
numClust=3;
dist_h='euclidean';
link='weighted';
hidx=clusterdata(PCA_TBTdata,'maxclust',numClust,'distance',dist_h,'linkage',link); % hierarchical clustering
%--------------------------------------------------------------------------
figure(11);
set(gcf,'Units','centimeter','Position',[5 0 18 20]);

alpha=0.05;
distribution='norm';
[~,Pos1]=P6_ConfidenceRegionFunction(PCA_TBTdata(find(hidx==1),1:3),alpha,distribution,[0.30,0.75,0.93],[0 0 0],'o',20);
hold on
[~,Pos2]=P6_ConfidenceRegionFunction(PCA_TBTdata(find(hidx==2),1:3),alpha,distribution,[0.30,0.75,0.93],[0 0 0],'o',20);
hold on
[~,Pos3]=P6_ConfidenceRegionFunction(PCA_TBTdata(find(hidx==3),1:3),alpha,distribution,[0.30,0.75,0.93],[0 0 0],'o',20);
hold on
  for m=1:3
      clear pos
      switch m
          case 1
              pos=Pos1;
          case 2
              pos=Pos2;
          case 3
              pos=Pos3;
      end
      for n=1:size(pos,1)
          id2=SecondSelectedFNBID(find(hidx==m),1);
          text(pos(n,1),pos(n,2),pos(n,3),[' FNB' num2str(id2(n,1))]);
      end
  end
%%
str1=sprintf('PC1(%4.2f%%)',explained(1,1));
str2=sprintf('PC2(%4.2f%%)',explained(2,1));
str3=sprintf('PC3(%4.2f%%)',explained(3,1));
xlabel(str1)
ylabel(str2)
zlabel(str3)
view(39,12)
set(gca,'Fontname', 'Arial','Fontsize', 10,'FontWeight','bold');

disp(0)
%% ====================================================================================
figure(12);
set(gcf,'Units','centimeter','Position',[5 0 18 20]);
scatter(0,0,0.01,'filled','o','MarkerFaceColor',[0.30,0.75,0.93],'MarkerEdgeColor',[0 0 0]) %突变体中非耐盐

text(0,5,'First Screening salt-tolerant FNB number:')
h=0;
count=0;
for m=1:size(FirstSelectedFNBID,1)
  text(0+0.55*count,4.5-h,num2str(FirstSelectedFNBID(m,1)))
  count=count+1;
  if count==15
      h=h+0.5; 
      count=0;
  end
  
end
text(0,2.5,'Second Screening salt-tolerant FNB number:')
h=0;
count=0;
for m=1:size(SecondSelectedFNBID,1)
  text(0+0.55*count,2-h,num2str(SecondSelectedFNBID(m,1)))
  count=count+1;
  if count==15
      h=h+0.5; 
      count=0;
  end
  
end

set(gca,'yLim',[0,5.5],'xTick',xtick,'XLim',[-0.5,8],'yTick',ytick,'Fontname', 'Arial',...
    'Fontsize', 10,'FontWeight','bold');
end

%% =========================================================================
function [TBT_Class_Result,wt,id,Bb]=Class_TBT_Select(Data)
%----------------------------------------------------
   D3_data=Data;
   [~,nMetaCount]=size(D3_data);
   nSample=36;
   TrainData=nan(nSample,nMetaCount);
   species=strings(nSample,1);
   
   
   GYH_Bio = normlization(D3_data,2);  %标准化
   TrainData2=GYH_Bio;
    TrainData=TrainData2(size(D3_data,1)-nSample+1:end,:);
   for m=1:nSample
       if m<=18
           species(m,1)='CK';
       else
           species(m,1)='Stress';
       end
   end
   
   sta_tab = tabulate(species);
   typenum = size(sta_tab,1);
   nsp = length(species);
   allabel = zeros(typenum,1);
   for i = 1:typenum
       currname = sta_tab(i,1);
       allabel(ismember(species,currname)) = i;
   end
   train_rate = 0.6; %
   train_num = nsp*train_rate;
   rand('seed',1); %随机种子
   random_sample=randperm(nsp);
   traindata = TrainData(random_sample(1:train_num),:);
   traingt = allabel(random_sample(1:train_num));
   testdata = TrainData(random_sample(train_num+1:end),:);
   testgt = allabel(random_sample(train_num+1:end));
     
   TBT_Data=GYH_Bio(1:end-nSample,:);  

    [kappa_cf, OA_cf, AA_cf,CA,curr_Allset_hat,curr_train_hat,wt,id,Bb]=SC_Classify_2022(TBT_Data,traindata,traingt,testdata,testgt);
   
     TBT_Class_Result=curr_Allset_hat(:,1);

     clear curr_Allset_hat;
%    end
end

%% ========================================================================
function data = normlization(data, choose)
% 数据归一化
  if choose==0
      % 不归一化
      data = data;
  elseif choose==1
      % Z-score归一化
      data = bsxfun(@minus, data, mean(data));
      data = bsxfun(@rdivide, data, std(data));
  elseif choose==2
      % 最大-最小归一化处理
      [data_num,~]=size(data);
      data=(data-ones(data_num,1)*min(data))./(ones(data_num,1)*(max(data)-min(data)));
  end
end
%% ========================================================================
function [kappa_cf, OA_cf, AA_cf, everyAA, Allsamples_hat,train_hat,weights,iranked,Bd] = SC_Classify_2022(alldata,traindata,traingt,testdata,testgt)

    traindata = double(traindata);
    testdata = double(testdata);
    alldata = double(alldata);
    traingt = uint8(traingt);
    testgt = uint8(testgt);

    weights=[];
    iranked=[];
    Bd=[];

    %%
    nTree = 100; 
    adnum=10;%
    traingt = traingt + adnum;  %
    RF_B = TreeBagger(nTree,traindata,traingt,'oobvarimp','on','MinLeafSize',1);  %增加特征重要性  %MinLeafSize
    test_hat = predict(RF_B, testdata);
    test_hat = str2num(cell2mat(test_hat));
    test_hat = uint8(test_hat - adnum);

    Allsamples_hat = predict(RF_B, alldata);
    Allsamples_hat = str2num(cell2mat(Allsamples_hat));
    Allsamples_hat = uint8(Allsamples_hat - adnum);

    train_hat = predict(RF_B, traindata);
    train_hat = str2num(cell2mat(train_hat)); %%
    train_hat = uint8(train_hat - adnum);  %%
    % 特征重要性
    weights=RF_B.OOBPermutedVarDeltaError;
    [Bd,iranked] = sort(weights,'descend'); % B:为各个特征的重要性得分（降序排列）;
    % iranked: 按照重要性的波段排序 weights 为各波段对应的权重
    [tt_CF_M, ~] = confusionmat(testgt,test_hat); %
    [kappa_cf, OA_cf, AA_cf,everyAA] = Output_Kappa(tt_CF_M); %
end
%% ========================================================================
function [h,PosT]=P6_ConfidenceRegionFunction(xdat,alpha,distribution,color1,color2,Marker,MarkerSize)
%   绘制置信区域（区间、椭圆区域或椭球区域）
%   CopyRight：xiezhh（谢中华）
%   https://www.ilovematlab.cn/thread-122092-1-1.html?_dsign=73a10ee6
% 设定参数默认值
if nargin == 1
    distribution = 'norm';
    alpha = 0.05;
elseif nargin == 2
    if ischar(alpha)
        distribution = alpha;
        alpha = 0.05;
    else
        distribution = 'norm';
    end
end
% 判断参数取值是否合适
if ~isscalar(alpha) || alpha>=1 || alpha<=0
    error('alpha的取值介于0和1之间')
end
if ~strncmpi(distribution,'norm',3) && ~strncmpi(distribution,'experience',3)
    error('分布类型只能是正态分布（''norm''）或经验分布（''experience''）')
end
% 检查数据维数是否正确
[m,n] = size(xdat);
p = min(m,n);  % 维数
if ~ismember(p,[1,2,3])
    error('应输入一维、二维或三维样本数据,并且样本容量应大于3')
end
% 把样本观测值矩阵转置，使得行对应观测，列对应变量
if m < n
    xdat = xdat';
end
xm = mean(xdat); % 均值
n = max(m,n);  % 观测组数
% 分情况绘制置信区域
switch p
    case 2    % 二维情形（置信椭圆）
        x = xdat(:,1);
        y = xdat(:,2);
        s = inv(cov(xdat));  % 协方差矩阵
        xd = xdat-repmat(xm,[n,1]);
        rd = sum(xd*s.*xd,2);
        if strncmpi(distribution,'norm',3)
            r = chi2inv(1-alpha,p);
        else
            r = prctile(rd,100*(1-alpha));

        end
        scatter(x(rd<=r),y(rd<=r),MarkerSize,Marker,'filled','MarkerEdgeColor',color2,'MarkerFaceColor',color1)  % 画样本散点

        hold on
        scatter(x(rd>=r),y(rd>=r),MarkerSize,Marker,'filled','MarkerEdgeColor',color2,'MarkerFaceColor',color1)  % 画样本散点

        h = ellipsefig(xm,s,r,1,color1*0.95);  % 绘制置信椭圆*0.25
        hold off;
    case 3    % 三维情形（置信椭球）
        x = xdat(:,1);
        y = xdat(:,2);
        z = xdat(:,3);
        s = inv(cov(xdat));  % 协方差矩阵
        xd = xdat-repmat(xm,[n,1]);
        rd = sum(xd*s.*xd,2);
        if strncmpi(distribution,'norm',3)
            r = chi2inv(1-alpha,p);
           % TitleText = '置信椭球（基于正态分布）';
        else
            r = prctile(rd,100*(1-alpha));
           % TitleText = '置信椭球（基于经验分布）';
        end
        plot3(x(rd<=r),y(rd<=r),z(rd<=r),'.','MarkerSize',16)  % 画样本散点
        hold on
        plot3(x(rd>r),y(rd>r),z(rd>r),'r+','MarkerSize',10)  % 画样本散点
        plot3(xm(1),xm(2),xm(3),'k+');  % 椭球中心
        h = ellipsefig(xm,s,r,2,color1*0.95);  % 绘制置信椭球
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        %title(TitleText)
        hidden off;
        hold off;
        PosT=[x,y,z];
end
end
%--------------------------------------------------
%   子函数：用来绘制置信区域（椭圆或椭球）
%--------------------------------------------------
function  h = ellipsefig(xc,P,r,tag,color)
        [V, D] = eig(P);
        if tag == 1
            aa = sqrt(r/D(1));
            bb = sqrt(r/D(4));
            t = linspace(0, 2*pi, 60);
            xy = V*[aa*cos(t);bb*sin(t)];  % 坐标旋转
            h = plot(xy(1,:)+xc(1),xy(2,:)+xc(2),'--', 'Color',color*0.8, 'linewidth', 1.5);
        else
            aa = sqrt(r/D(1,1));
            bb = sqrt(r/D(2,2));
            cc = sqrt(r/D(3,3));
            [u,v] = meshgrid(linspace(-pi,pi,30),linspace(0,2*pi,30));
            x = aa*cos(u).*cos(v);
            y = bb*cos(u).*sin(v);
            z = cc*sin(u);
            xyz = V*[x(:)';y(:)';z(:)'];  % 坐标旋转
            x = reshape(xyz(1,:),size(x))+xc(1);
            y = reshape(xyz(2,:),size(y))+xc(2);
            z = reshape(xyz(3,:),size(z))+xc(3);
            h = mesh(x,y,z);  % 绘制椭球面网格图
        end
    end
%%
function [k, po, aa,everyAA] = Output_Kappa(varargin)
% 原始的Kappa函数见kappa.m文件，原始文件是以混淆矩阵为输入，直接在command Window打印出多个结果，函数无返回值
% 原始文件的官网链接：http://cn.mathworks.com/matlabcentral/fileexchange/15365-cohen-s-kappa/content/kappa.m
% 为了能得到Kappa返回值，进行了函数的改写，如果想要得到更多的返回值，参考原始文件再进行输出改写
% 返回值k为混淆矩阵对应的预测值与真实值的Kappa系数；po为混淆矩阵对应的预测值的总体精度；aa为平均分类精度


args=cell(varargin);
nu=numel(args);
if isempty(nu)
    error('Warning: Matrix of data is missed...')
elseif nu>3
    error('Warning: Max three input data are required')
end
default.values = {[],0,0.05};
default.values(1:nu) = args;
[x w alpha] = deal(default.values{:});

if isempty(x)
    error('Warning: X matrix is empty...')
end
if isvector(x)
    error('Warning: X must be a matrix not a vector')
end
if ~all(isfinite(x(:))) || ~all(isnumeric(x(:)))
    error('Warning: all X values must be numeric and finite')
end   
if ~isequal(x(:),round(x(:)))
    error('Warning: X data matrix values must be whole numbers')
end

m=size(x); % m为行列数 m(1)为行数，m(2)为列数
if ~isequal(m(1),m(2))
    error('Input matrix must be a square matrix')
end
if nu>1 %eventually check weight
    if ~isscalar(w) || ~isfinite(w) || ~isnumeric(w) || isempty(w)
        error('Warning: it is required a scalar, numeric and finite Weight value.')
    end
    a=-1:1:2;
    if isempty(a(a==w))%check if w is -1 0 1 or 2
        error('Warning: Weight must be -1 0 1 or 2.')
    end
end
if nu>2 %eventually check alpha
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
end
clear args default nu

m(2)=[]; % 将m(2)删去，即得到m变为原来的m的行数，亦即m为混淆矩阵的行数
f=diag(ones(1,m)); %unweighted
n=sum(x(:)); %Sum of Matrix elements
x=x./n; %proportion
r=sum(x,2); %rows sum
s=sum(x); %columns sum
Ex=r*s; %expected proportion for random agree
pom=sum(min([r';s]));
po=sum(sum(x.*f));
pe=sum(sum(Ex.*f));
k=(po-pe)/(1-pe);
km=(pom-pe)/(1-pe); %maximum possible kappa, given the observed marginal frequencies
ratio=k/km; %observed as proportion of maximum possible
sek=sqrt((po*(1-po))/(n*(1-pe)^2)); %kappa standard error for confidence interval
ci=k+([-1 1].*(abs(-realsqrt(2)*erfcinv(alpha))*sek)); %k confidence interval
wbari=r'*f;
wbarj=s*f;
wbar=repmat(wbari',1,m)+repmat(wbarj,m,1);
a=Ex.*((f-wbar).^2);
var=(sum(a(:))-pe^2)/(n*((1-pe)^2)); %variance
z=k/sqrt(var); %normalized kappa
p=(1-0.5*erfc(-abs(z)/realsqrt(2)))*2;

% 计算平均分类精度AA(Average accuracy)：各个类别被正确分类的百分比的均值
% 通过混淆矩阵计算：主对角线上的值与对应的行和的比值的均值，即为AA
for i = 1 : m
    everyAA(i) = x(i,i) / r(i);
end
aa = mean(everyAA);
end

