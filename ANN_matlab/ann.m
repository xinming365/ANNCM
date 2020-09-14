TotalDataBaseRead = xlsread('数据文件.xls','B4:G19');
X1 = TotalDataBaseRead(:,1);
X2 = TotalDataBaseRead(:,2);
X3 = TotalDataBaseRead(:,3);
X4 = TotalDataBaseRead(:,4);
Y1 = TotalDataBaseRead(:,5);
Y2 = TotalDataBaseRead(:,6);
 
InputsAll = [X1,X2,X3,X4];
OutputOne = Y1;
OutputTwo = Y2;
Groups = combnk(1:8,3)'; % 产生不同排列组合，搜索寻找不同参数之间的关系，并将根据相关性程度的不同进行参数分析的筛选
[rowofGroups,columnofGroups]=size(Groups);
 
for iii = 1:columnofGroups
    InputsIndex = Groups(:,iii);
    Inputs = InputsAll(:,InputsIndex);
    TotalDatabaseO = [Inputs OutputOne]; 
    [row,column]=size(TotalDatabaseO);
    TrainingNumber=round(row*4/5);
    NNHLMax=fix((TrainingNumber-1)/(column+1));% 确定隐层网络的最大神经元个数
    NDistribu=80; % 确定数据混排次数
    TotalDataBaseN=zeros(row,column,NDistribu,NNHLMax); 
    mrand1=zeros(NDistribu,NNHLMax); brand1=zeros(NDistribu,NNHLMax); rrand1=zeros(NDistribu,NNHLMax);
    mrand2=zeros(NDistribu,NNHLMax); brand2=zeros(NDistribu,NNHLMax); rrand2=zeros(NDistribu,NNHLMax);
    OutCell=cell(17,13);
    EvaluateFun=zeros(NDistribu,NNHLMax); MinEvaluateFunoftime=zeros(1,NNHLMax); indexoftime=zeros(1,NNHLMax);
    index=zeros(row,1,NDistribu,NNHLMax);
    
    for k = 1:NNHLMax
        for kk = 1:NDistribu
            if k==1
                AddNumber=rand(row,1);
                [AddNumberNew,index(:,1,kk,k)]=sort(AddNumber);
                TotalDataBaseN(:,:,kk,k)=TotalDatabaseO(index(:,1,kk,k),:);
            end
            TotalDataBaseN(:,:,kk,k)=TotalDataBaseN(:,:,kk,1);
            index(:,:,kk,k)=index(:,:,kk,1);
            PropertiesInNew=TotalDataBaseN(:,(1:(column-1)),kk,k)';
            PropertiesOutNew=TotalDataBaseN(:,column,kk,k)';
%             [PropertiesInNewn,PropertiesInNews]= mapminmax(PropertiesInNew);
%             [PropertiesOutNewn,PropertiesOutNews]=mapminmax(PropertiesOutNew); 
            iitst=3:5:row;
            iitr=[1:5:row 2:5:row 4:5:row 5:5:row];
            test.P=PropertiesInNew(:,iitst);test.T=PropertiesOutNew(:,iitst);
            ptr=PropertiesInNew(:,iitr);ttr=PropertiesOutNew(:,iitr);
            net=newff(ptr,ttr,k,{},'trainbr');
            net.trainParam.epochs = 800;
            net.trainParam.goal = 1e-8;

            net=init(net);
            [net,tr]=train(net,ptr,ttr);
            PropertiesOutSimulateTr=sim(net,ptr);
            PropertiesOutSimulatets=sim(net,test.P);
%             PropertiesOutSimulateTr=mapminmax('reverse',PropertiesOutSimulatenTr,PropertiesOutNews);
%             PropertiesOutSimulatets=mapminmax('reverse',PropertiesOutSimulatents,PropertiesOutNews);
            [mrand1(kk,k),brand1(kk,k),rrand1(kk,k)]=postregnopic(PropertiesOutSimulateTr,PropertiesOutNew(:,iitr)); % Éñ¾­ÍøÂç¹¹½¨¹ý³ÌÖÐÆÀ¼ÛÆä²¶×½ÄÜÁ¦
            [mrand2(kk,k),brand2(kk,k),rrand2(kk,k)]=postregnopic(PropertiesOutSimulatets,PropertiesOutNew(:,iitst));
            if isnan(rrand1(kk,k))
                rrand1(kk,k)=0;
            elseif rrand1(kk,k)==Inf
                rrand1(kk,k)=0;
            elseif rrand1(kk,k)==-Inf
                rrand1(kk,k)=0;
            end
            if isnan(rrand2(kk,k))
                rrand2(kk,k)=0;
            elseif rrand2(kk,k)==Inf
                rrand2(kk,k)=0;
            elseif rrand2(kk,k)==-Inf
                rrand2(kk,k)=0;
            end
            EvaluateFun(kk,k)=(((abs(mrand1(kk,k)-1)+(1-rrand1(kk,k)))+(abs(mrand2(kk,k)-1)+(1-rrand2(kk,k)))))+(abs((abs(mrand1(kk,k)-1)+(1-rrand1(kk,k)))-(abs(mrand2(kk,k)-1)+(1-rrand2(kk,k)))));
        end
        [MinEvaluateFunoftime(k),indexoftime(k)]=min(EvaluateFun(:,k));
    end
    [MinEvaluateFunofNeurons,indexofneuron]=min(MinEvaluateFunoftime(:));
    TotalDataBaseNFinal=TotalDataBaseN(:,:,indexoftime(indexofneuron),indexofneuron);
    PropertiesInNewFinal=TotalDataBaseNFinal(:,(1:(column-1)))';
    PropertiesOutNewFinal=TotalDataBaseNFinal(:,column)';
%     [PropertiesInNewFinaln,PropertiesInNewFinals]=mapminmax(PropertiesInNewFinal);
%     [PropertiesOutNewFinaln,PropertiesOutNewFinals]=mapminmax(PropertiesOutNewFinal);
    iitstFinal=3:5:row;
    iitrFinal=[1:5:row 2:5:row 4:5:row 5:5:row];
    testFinal.P=PropertiesInNewFinal(:,iitstFinal);testFinal.T=PropertiesOutNewFinal(:,iitstFinal);
    ptrFinal=PropertiesInNewFinal(:,iitrFinal);ttrFinal=PropertiesOutNewFinal(:,iitrFinal);
    netFinal=newff(ptrFinal,ttrFinal,indexofneuron,{},'trainbr');
    netFinal.trainParam.epochs = 800;
    netFinal.trainParam.goal = 1e-8;
%     netFinal=init(netFinal);
    [netFinal,trFinal]=train(netFinal,ptrFinal,ttrFinal);
    PropertiesOutNewFinalSimulateTr=sim(netFinal,ptrFinal);
    PropertiesOutNewFinalSimulatets=sim(netFinal,testFinal.P);
%     PropertiesOutNewFinalSimulateTr=mapminmax('reverse',PropertiesOutNewFinalSimulatenTr,PropertiesOutNewFinals);
%     PropertiesOutNewFinalSimulatets=mapminmax('reverse',PropertiesOutNewFinalSimulatents,PropertiesOutNewFinals);
    PropertiesOutNewFinalTr=PropertiesOutNewFinal(:,iitrFinal);
    PropertiesOutNewFinalts=PropertiesOutNewFinal(:,iitstFinal);
    Erroroftesting=PropertiesOutNewFinalts-PropertiesOutNewFinalSimulatets;
    ErroroftestingPercentage=Erroroftesting./PropertiesOutNewFinalts;
    FinalProperty=[PropertiesOutNewFinalTr PropertiesOutNewFinalts];
    FinalSimulate=[PropertiesOutNewFinalSimulateTr PropertiesOutNewFinalSimulatets];
    Errorofwhole=FinalProperty-FinalSimulate;
    ErrorofwholePercentage=Errorofwhole./FinalProperty;
    
    % Testing
    MeanErroroftesting=mean(abs(Erroroftesting));
    StdvErroroftesting=std(abs(Erroroftesting));
    MeanErroroftestingPercentage=mean(abs(ErroroftestingPercentage));
    StdvErroroftestingPercentage=std(abs(ErroroftestingPercentage));
    
    % Whole
    MeanErrorofwhole=mean(abs(Errorofwhole));
    StdvErrorofwhole=std(abs(Errorofwhole));
    MeanErrorofwholePercentage=mean(abs(ErrorofwholePercentage));
    StdvErrorofwholePercentage=std(abs(ErrorofwholePercentage)); 
    
    strPI = cell(rowofGroups,1);
    for jjj=1:rowofGroups
        switch Groups(jjj,iii)
            case 1
                strPI(jjj) = {'X1'};
            case 2
                strPI(jjj) = {'X2'};
            case 3
                strPI(jjj) = {'X3'};
            case 4
                strPI(jjj) = {'X4'};
        end
    end
    strPO = 'Y1'; % Can change to Y2
    [m3,b3,r3]=postregsfinal(PropertiesOutNewFinalSimulateTr,PropertiesOutNewFinalSimulatets,PropertiesOutNewFinalTr,PropertiesOutNewFinalts);
    FigureCrit=abs(1-m3)+(1-r3);
    nameoffileTr=['(',num2str(FigureCrit),')',' ','Predict',' ',strPO,' ','from',' ',char(strPI(1)),'(training)','.fig']; % Change name
    nameoffileTs=['(',num2str(FigureCrit),')',' ','Predict',' ',strPO,' ','from',' ',char(strPI(1)),'(testing)','.fig']; % Change name
    nameoffileWhole=['(',num2str(FigureCrit),')',' ','Predict',' ',strPO,' ','from',' ',char(strPI(1)),'(whole set)','.fig']; % Change name        
    nameoffile2=['(',num2str(FigureCrit),')',' ','Predict',' ',strPO,' ','from',' ',char(strPI(1)),'(whole set)','.xls'];
    h3=gcf;
    saveas(h3,nameoffileWhole,'fig')
    delete(h3)
    [m1,b1,r1]=postreglabel(PropertiesOutNewFinalSimulateTr,PropertiesOutNewFinalTr); % Label the axis in detail
    h1=gcf;
    saveas(h1,nameoffileTr,'fig')
    delete(h1)
    [m2,b2,r2]=postreglabel(PropertiesOutNewFinalSimulatets,PropertiesOutNewFinalts); % Label the axis in detail
    h2=gcf;
    saveas(h2,nameoffileTs,'fig')
    delete(h2) 
    OutCell(1,:)={'Predicted Error of Testing','Percentage Error of Testing','Predicted Error of Whole Set','Percentage Error of Whole Set','Mean Error of Testing (modul)','Error Standard Deviation of Testing (modul)','Mean Percentage Error of Testing (modul)','Percentage Error Standard Deviation of Testing (modul)','Mean Error of Whole (modul)','Error Standard Deviation of Whole (modul)','Mean Percentage Error of Whole (modul)','Percentage Error Standard Deviation of Whole (modul)','Number of Neurons in Hidden Layer'};
    OutCell(2:4,1)=num2cell(Erroroftesting)';
    OutCell(2:4,2)=num2cell(ErroroftestingPercentage)';
    OutCell(2:17,3)=num2cell(Errorofwhole)';
    OutCell(2:17,4)=num2cell(ErrorofwholePercentage)';
    OutCell(2,5:13)={MeanErroroftesting,StdvErroroftesting,MeanErroroftestingPercentage,StdvErroroftestingPercentage,MeanErrorofwhole,StdvErrorofwhole,MeanErrorofwholePercentage,StdvErrorofwholePercentage,indexofneuron};
    xlswrites(nameoffile2,OutCell,'A1:M18')    
end
