#include "head.h"
#include "source.cpp"
#include "mat.cpp"

/**
 * 可以多次测试得到最好的p,q 
 */
int Calculate_pq(vector<Double> data)
{

    freopen("r_arma_plus.xls", "w", stdout);

    Double mean; //输入数据的均值
    vector<Double> AutoCor;//自相关系数AutoCorrelation
    vector<Double> BiasCor;//偏相关系数

    AutoCor = getAutoCor(data); //得到的自相关系数
	BiasCor = getBiasCor(AutoCor); // 得到偏相关系数
	//输出到excel文件中，直接用图表显示，便于得到p,q
	for(int k=0;k<AutoCor.size();k++){
        cout<<AutoCor[k]<<"\t";
    }
	cout<<endl;
    for(int k=0;k<BiasCor.size();k++){
        cout<<BiasCor[k]<<"\t";
    }

    return 0;
}



/**
 *根据p得到数据矩阵形式
 */
vector<vector<Double> > y; //y = [data[p+1, ..., n]]
vector<vector<Double> > x;
/**
  [ data[p, ..., 1] ]
  [ data[p+1, ..., 2] ]  
  [ data[p+2, ..., 3] ]
  .
  .
  .
  [ data[n-1, ..., n-p] ]
 */


void formatData(vector<Double> data,int p){
	vector<Double> tmpy;
	for(int i=p;i<data.size();i++){
		tmpy.push_back(data[i]);
		vector<Double> tmp;
		for(int j=i-1;j>=i-p;j--){
			tmp.push_back(data[j]);
		}
		x.push_back(tmp);
	}
	y.push_back(tmpy);
	y = t(y);
}

/**
 *最小二乘法求参数
 * a = inv(t(x) _*_ x) _*_ t(x) _*_ Y
 * e = sum(a) / (n-p)
 */
vector<Double> LeastSquares(vector<Double> data,int p){
	formatData(data,p);	
	
	vector<vector<Double> > a, tx,invx,tmp;
	tx = t(x);
	invx = inv(mulMat(tx, x));
	a = mulMat(mulMat(invx,tx), y);
	a = t(a);
	return a[0];
}


/*
 *3、由上述估计的 AR(PP)模型递推计算残差列 { ε[t] }
 */
vector<Double> getBiasSeries(vector<Double> data,vector<Double> a,int p){
	vector<Double> calPN(data.begin(),data.begin()+p);
	
	for(int i=p;i<data.size();i++){
		Double s = 0;
		int t = calPN.size();
		for(int j=0;j<p;j++){
			s += a[j] * calPN[t-j-1];
		}
		calPN.push_back(s);
	}
	
//	cout<<calPN.size()<<endl;
	vector<Double> var;
	//计算残差
	for(int i=p;i<calPN.size();i++){
		var.push_back(data[i] - calPN[i]);
	}
	/**
	cout<<"残差"<<endl;
	for(int k=0;k<var.size();k++){
        cout<<var[k]<<"\t";
    }
    **/
    return var;

}



/**
 *自回归逼近方法求参数
 *1、首先利用原始数据 x1 , ... , xn 作高阶自回归滑动模型 AR(pp) 的拟合
 *	用较高阶的 AR(PP) 模型 取 PP >> p,q才能较好的拟合原始数据序列
 *2、确定pp后，估计AR(pp)中参 数α1，α2， ... , αn使其满足x[t] = α[1]*x[t−1] + α[2]x[t−2] + ... + α[pp] x[t−pp] + ε[t]
 *3、由上述估计的 AR(PP)模型递推计算残差列 { ε[t] }
 *4、再视残差列 { ε[t] } 为独立序列 利用线性回归模型
 *5、利用矩阵可得 α 与 β 的最小二乘估计
 *
 * 这里选择 p = 7， q = 20, pp = 35，原则应该选择q=20，考虑到计算矩阵是可能会超出C++语言的精度，如果超出，可以选择q = 8进行计算
 */


/**
X = [ x[pp+q], x[pp+q-1], ..., x[pp+q-p+1] ]
	[ x[pp+q+1], x[pp+q], ..., x[pp+q-p+2] ]
	.
	.
	.
	.
	.
	[ x[n-1], x[n-2], ..., x[n-p] ]


E = [ e[p+q], e[p+q-1], ..., e[p+1] ]
	[ e[p+q+1], e[p+q], ..., e[p+2] ]
	.
	.
	.
	.
	.
	[ e[n-1],e[n-2], ..., e[n-q] ]


**/
vector<vector<Double> > X,E,xs;
void formatData(vector<Double> data,vector<Double> bias,int p, int q,int pp){
	int n = data.size();
	for(int i=pp;i<n-q;i++){
		vector<Double> tmp;
		for(int j=0;j<p;j++){
			tmp.push_back(data[i+q-j]);
		}
		X.push_back(tmp);
	}
	
	for(int i=0;i<bias.size()-q;i++){
		vector<Double> tmp;
		for(int j=0;j<q;j++){
			tmp.push_back(bias[i+q-j]);
		}
		E.push_back(tmp);
	}
	
	vector<Double> tmpx;
	for(int i=pp+q;i<n;i++){
		tmpx.push_back(data[i]);
	}
	
	xs.push_back(tmpx);
	xs = t(xs);
}

/**
 *得到参数,a,b
 */

vector<Double> getParm_ab(vector<Double> data,vector<Double> bias,int p, int q,int pp){
	formatData(data,bias,p,q,pp);	
	
	vector<vector<Double> > ab, tx,invxe,tmp, r_XE, c_XE;
	c_XE = ConCols(t(X),t(E));
	r_XE = ConRows(X,E);
	
	invxe = inv(mulMat(c_XE,r_XE));
	tmp = mulMat(invxe,c_XE);
	
	ab = mulMat(tmp,xs);

	return t(ab)[0];	
}


/**
 *x[t] = sum[j: 0...p]{a[j]*data[t-j]} + e
 *
 *
 *ARMA ( p , q ) 模型的拟合检验
 *检验 ARMA ( p , q ) 模型与检验 AR ( p ) 模型 MA (q ) 模型的方法基本相同 都是检验其拟合残差序列是否为独立序列 不同的是 各自获得拟合残
 *差序列序列时使用各自不同的拟合模型而已
 *即实际上
 *检验 x[t] 是否为
 *ARMA ( p , q ) 时 只需检验残差列 { ε t } 是否独立序列即可 而残差列的估计
 *值 { ε ˆ t } 可由样本值 x1 , ... , xn 计算得出
 *判断 { ε ˆ t } 是否独立序列
 *最后再利用判别独立序列的方法
 *则认为 { x t } 为 ARMA ( p , q ) 序列
 *若是,否则认 { x t } 不是 ARMA ( p , q ) 序列
 *
 *1、提出模型结社H0
 *2、带入得到的参数，得到x t = ∑[i,..,p]α[i]*x[t−i] + ε[t] − ∑[k,...,q]β[k]*ε[t−k]
 *3、由上式计算残差
 *4、由 { ε[t] } 求自协方差函数
 *检验原则：
 *若 { ρ(ε), k = 0,1,2, ... , n } 中 约 有 68.3% 的 点 落 在 纵 坐 标ρ = ± 1 / n 内
 *约 有 95.4% 的 点 落 在 纵 坐 标 ρ = ± 2 / n 内
 *( ε[p+1] , ε[p+2] ,..., ε[n]) 为独立序列样本值 此时接受H0 否则拒绝H0
 */

int calPQ_N(vector<Double> data,vector<Double> data_var,vector<Double> a,vector<Double> b,int p,int q){
	
	int n = data.size();	
	//得到残差第二项x[t] − ∑α[k]*x[t−k] t = [p, ... ,N]
	vector<Double> calPN(data.begin(),data.begin()+p);
		
	Double res = 0;
	for(int i=p;i<n;i++){
		Double s = 0;
		int t1 = data.size();
		for(int j=0;j<p;j++){
			s += a[j] * data[t1-j-1];
		}
		int t2 = data_var.size();
		for(int j=0;j<q;j++){
			s -= b[j] * data_var[t2-j-1];
		}		
		calPN.push_back(s);
	}	
	
	//cout<<calPN.size()<<endl;
	vector<Double> var;
	//计算残差
	for(int i=p;i<calPN.size();i++){
		var.push_back(data[i] - calPN[i]);
	}
	
	vector<Double> varpq;//得到ε[p,.., n]的残差，求自协方差，判断H0是否独立
	for(int t=p;t<n;t++){
		Double tmp = 0;
		for(int j=0;j<q;j++){
			tmp += b[j]*var[t-j-1];
		}
		var.push_back(tmp - var[t-p]);
		varpq.push_back(tmp - var[t-p]);
	}
	
	vector<Double> Cor = getAutoCor(varpq);
	
	freopen("rvar_ARMA_p7_q_8.xls", "w", stdout);
	cout<<"自相关系数:"<<endl;
	for(int k=0;k<Cor.size();k++){
        cout<<Cor[k]<<"\t";
    }
    cout<<endl;
	
	//检验是否有 68.3% 的 点 落 在 纵 坐 标ρ = ± 1 / n 内
    //约 有 95.4% 的 点 落 在 纵 坐 标 ρ = ± 2 / n 内
    int k1 = 0,k2 = 0;
    Double p1 = 1.0 / Cor.size(),p2 = 2.0 / Cor.size();
    for(int k=0;k<Cor.size();k++){
    	if(Cor[k] >= -p1 && Cor[k] <= p1) k1++;	
    	if(Cor[k] >= -p2 && Cor[k] <= p2) k2++;	
    }
    cout<<"ρ = ± 1 / n："<<k1*1.0 / Cor.size()<<endl;
    cout<<"ρ = ± 2 / n："<<k2*1.0 / Cor.size()<<endl;
	return 0;	
	
}


/**
 *根据模型预测以后的数据，k表示第k个数据，这里k大于n,注意时间序列，必须先预测得到n，才能得到n+1，
 *如果给出的k>n，会预测[n,k]的所有位置，并添加大原数据上。
 *vector<Double> &data : 原始数据
 *vector<Double> &data_var: 根据数据得到的残差
 *a、b、e : 模型参数
 *k : 求的第k个数据
 */

Double predict(vector<Double> &data,vector<Double> &data_var,vector<Double> a,vector<Double> b,int p,int q,int k){
	Double res = 0;
	for(int i=data.size();i<k;i++){
		Double s = 0;
		int t1 = data.size();
		for(int j=0;j<p;j++){
			s += a[j] * data[t1-j-1];
		}
		int t2 = data_var.size();
		for(int j=0;j<q;j++){
			s -= b[j] * data_var[t2-j-1];
		}		
		data.push_back(s);
	}
	return data[k-1];
}




Double xx[] = {871.5, 897.1, 904.3, 919.2, 935.0, 950.0, 965.0, 981.0,1028.0,1047.0,1061.0,1075.0,
              1086.0,1094.0,1110.0,1112.0,1125.0,1151.1,1159.4,1180.0,1195.6,1227.2,1243.6,
              1256.0,1128.0,1292.0,1296.0,1298.0,1302.1,1309.4,1317.0,1332.6,1235.2,1363.6,
              1385.1,1423.2,1456.4,1472.7,1488.0,1491.0,1501.0,1511.0,1520.0,1531.9,1538.6,1540.3
};

Double xd[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,18,19,20};

int main()
{
	vector<Double> data;
	int p=7,q=8,pp=20;//一定注意p，q的取值是通过数据计算后，估计出来的。
	//读入数据
	for(int i=0;i<46;i++){
		data.push_back(xx[i]);
	}
	//计算p,q,通过图像显示，选择，p = 7， q = 20, pp = 38
//	Calculate_pq(data);
	
	vector<Double> ta = LeastSquares(data,pp);
	cout<<"根据AR模型得到的参数ta个数:  "<<ta.size()<<endl;
	for(int i=0;i<ta.size();i++){
		cout<<"ta["<<i<<"] = "<<ta[i]<<endl;
	}
	
	//残差	
	vector<Double> bias = getBiasSeries(data,ta,pp);
	/**
	for(int i=0;i<bias.size();i++){
		cout<<"var["<<i<<"] = "<<bias[i]<<endl;
	}
	**/
	vector<Double> ab = getParm_ab(data,bias,p,q,pp);

	vector<Double> a(ab.begin(),ab.begin()+p);	
	vector<Double> b(ab.begin()+p,ab.begin()+p+q);
	cout<<"参数a个数:  "<<a.size()<<endl;
	for(int i=0;i<a.size();i++){
		cout<<"a["<<i<<"] = "<<a[i]<<endl;
	}
	cout<<"参数b个数:  "<<b.size()<<endl;
	for(int i=0;i<b.size();i++){
		cout<<"b["<<i<<"] = "<<b[i]<<endl;
	}
	
	calPQ_N(data,bias,a,b,p,q);
	
	cout<<predict(data,bias,a,b,p,q,47)<<endl;
	
	
	return 0;
}



















