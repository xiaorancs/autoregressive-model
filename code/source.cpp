#include "head.h"
/**
 *该文件中是两个算法都用到的函数
 */




/**
 *自协方差 AutoCov[k] = E((x[i] - u)(x[i-k] - u))
 *自相关系数 AutoCov[k] = AutoCov[k] / AutoCov[0]
 */
 
vector<Double> getAutoCov(vector<Double> data){
	//计算自相关系数矩阵
    int n = data.size();

    Double mean = 0; //数据的均值
    for(int i=0;i<n;i++){
        mean += data[i];
    }
    mean /= n;
	//cout<<"mean::"<<mean<<endl;
    //将每个数据都减去均值得到新的数据
    vector<Double> prodata;
    
    for(int i=0;i<n;i++){
        prodata.push_back(data[i] - mean);
		//cout<<"prodata[i] "<<prodata[i]<<endl;
    }

    vector<Double> AutoCov(n,0);//自协方差AutoCovariance
    for(int k=0;k<n;k++){
        for(int i=0;i<n-k;i++){
            AutoCov[k] += prodata[i] * prodata[i+k];
        }
        AutoCov[k] /= n - k;
    }

	return AutoCov;
}



vector<Double> getAutoCor(vector<Double> data){
    
	vector<Double> AutoCor,AutoCov;//自相关系数AutoCorrelation，注意下标从0开始
    AutoCov = getAutoCov(data);
    
    /**
    cout<<"AutoCor:自协方差：begin"<<endl;
    for(int k=0;k<AutoCov.size();k++){
    	cout<<AutoCov[k]<<" ";
    }
    cout<<"AutoCor:自协方差：end"<<endl;
    **/
    
    for(int k=0;k<data.size()-1;k++){
        AutoCor.push_back(AutoCov[k+1] / AutoCov[0]);
    }
	
	/**
	cout<<"AutoCor:自协方差系数：begin"<<endl;
    for(int k=0;k<AutoCor.size();k++){
    	cout<<AutoCor[k]<<" ";
    }
    cout<<"AutoCor:自协方差系数：end"<<endl;
    **/
    return AutoCor;
}


/**
 *得到偏相关系数BiasCor[k,k]
 *BiasCor[0,0] = AutoCor[0]
 *BiasCor[k,k] = (AutoCor[k-1] - sum[j:0...k-1]{AutoCor[k-j]*BiasCor[j,k-1]}) / (1 - sum[j:0...k-1]AutoCor[j]*BiasCor[j,k-1])
 *BiasCor[j,k] = BiasCor[j,k-1] - BiasCor[k,k]*BiasCor[k-j,k-1] j = 0...k
 *
 */
vector<Double> getBiasCor(vector<Double> AutoCor){
    //计算BiasCor[i,j],为了直接访问下标，首先初始化
    vector< vector<Double> > BiasCor;
    for(int i=0;i<AutoCor.size();i++){
        vector<Double> tmp(AutoCor.size(),0);
        BiasCor.push_back(tmp);
    }

    BiasCor[0][0] = AutoCor[0];

    for(int k=1;k<AutoCor.size();k++){
        BiasCor[k][k] = AutoCor[k];
        Double t1,t2;
        for(int j=0;j<=k-1;j++){
            t1 = AutoCor[k-j] * BiasCor[j][k-1];
            t2 = AutoCor[j] * BiasCor[j][k-1];

            BiasCor[j][k] = BiasCor[j][k-1] - BiasCor[k][k] * BiasCor[k-j][k-1];

        }
        BiasCor[k][k] = (BiasCor[k][k] - t1) / t2;
        for(int j=0;j<=k-1;j++){
            BiasCor[k][j] = BiasCor[j][k] = BiasCor[j][k-1] - BiasCor[k][k] * BiasCor[k-j][k-1];
        }
    }
    vector<Double> res;
    for(int k=0;k<AutoCor.size();k++){
        res.push_back(BiasCor[k][k]);
    }

    return res;

}
