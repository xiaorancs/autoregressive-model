#include "head.h"

/**
 * a = inv(t(x) _*_ x) _*_ t(x) _*_ Y
 * t(x) : 对x球转置
 * r(x) : 对x求逆矩
 * -*- : 矩阵乘法
 * inv(x): 矩阵的逆
 */

/**
 *矩阵转置
 */
vector<vector<Double> > t(vector<vector<Double> > x){
    //x的装置矩阵
    vector<vector<Double> > tx;
    //tx初始化便于直接访问下标,这是原矩阵的转置的形式
    for(int i=0;i<x[0].size();i++){
        vector<Double> tmp(x.size(),0);
        tx.push_back(tmp);
    }

    for(int i=0;i<x.size();i++){
        for(int j=0;j<x[0].size();j++){
            tx[j][i] = x[i][j];
        }
    }
    return tx;
}

/**
 *矩阵乘法
 */
vector<vector<Double> > mulMat(vector<vector<Double> > tx, vector<vector<Double> > x){
    vector<vector<Double> > res;
    //初始化结果矩阵的格式row(tx) X col(x)
    for(int i=0;i<tx.size();i++){
        vector<Double> tmp(x[0].size(),0);
        res.push_back(tmp);
    }

//    cout<<res.size()<<" "<<res[0].size()<<endl;
    for(int i=0;i<tx.size();i++){
        for(int j=0;j<x[0].size();j++){
            for(int k=0;k<x.size();k++){
               res[i][j] += tx[i][k] * x[k][j];
            }
        }
    }
    return res;
}

/**
 *矩阵的行列式，行列变化为上三角矩阵
 */

Double det(vector<vector<Double> > x){
	//只有一个元素
	//if(x.size() == 1 && x[0].size() == 1) return x[0][0];
	
	Double det = 1;
	//交换数组指定的两行，即进行行变换（具体为行交换）
    int iter = 0;  //记录行变换的次数（交换）
    for(int i=0;i<x.size();i++){
        if(x[i][i] == 0){
        	for(int j=i+1;j<x.size();j++){
	        	if(x[j][i] != 0){
	                swap(x[i],x[j]);//交换两行
    	            iter ++;
    	        }
    	    }
        }
        for(int k=i+1;k<x.size();k++){
            Double yin = -1 * x[k][i] / x[i][i] ;
            for(int u=0; u<x[0].size(); u++){
                x[k][u] = x[k][u] + x[i][u] * yin;
            }
        }
    }
	
	/**
   	cout<<"上三角矩阵："<<endl;
   	for(int i=0;i<x.size();i++){
        for(int j=0;j<x[0].size();j++){
            cout<<x[i][j]<<" ";
        }
        cout<<endl;
    }**/
	for(int i=0;i<x.size();i++){//求对角线的积 即 行列式的值
      det = det * x[i][i];
	}
	//行变换偶数次符号不变
	if(iter%2 == 1)  det= -det;

	return det;
}

/**
 *删除矩阵的第r行，第c列
 */
vector<vector<Double> > delMat(vector<vector<Double> > x,int r,int c){
	vector<vector<Double> > Ax;
	for(int i=0;i<x.size();i++){
		vector<Double> tmp;
		for(int j=0;j<x[0].size();j++){
			if(i != r && j != c) tmp.push_back(x[i][j]);			
		}
		if(i != r) Ax.push_back(tmp);
	}	
	return Ax;
}


/**
 *求矩阵的伴随矩阵
 */
vector<vector<Double> > A(vector<vector<Double> > x){
	vector<vector<Double> > tmp(x),res;
	
	//tx初始化便于直接访问下标,这是原矩阵的转置的形式
    for(int i=0;i<x.size();i++){
        vector<Double> tp(x[0].size(),0);
        res.push_back(tp);
    }
    
	for(int i=0;i<x.size();i++){
		for(int j=0;j<x[0].size();j++){
			tmp = x;
			tmp = delMat(tmp,i,j);
			res[i][j] = ((i+j)%2==0?1:-1) * det(tmp);
			
		}
	}
	return res;
}


/**
 *矩阵的逆
 */
vector<vector<Double> > inv(vector<vector<Double> > x){
	vector<vector<Double> > res = A(x);
	Double dets = det(x);
	for(int i=0;i<res.size();i++){
		for(int j=0;j<res[0].size();j++){
			res[i][j] /= dets;
		}
	}
	return res;
}


/**
 *合并两个行相同的矩阵
 */ 
vector<vector<Double> > ConRows(vector<vector<Double> > x, vector<vector<Double> > y){
	//行相同，添加列
	for(int i=0;i<y.size();i++){
		for(int j=0;j<y[0].size();j++){
			x[i].push_back(y[i][j]);
		}
	}
	return x;
}

/**
 *合并两个列相同的矩阵
 */ 
vector<vector<Double> > ConCols(vector<vector<Double> > x, vector<vector<Double> > y){
	//列相同，添加行
	for(int i=0;i<y.size();i++){
		vector<Double> row;
		for(int j=0;j<y[0].size();j++){
			row.push_back(y[i][j]);
		}
		x.push_back(row);
	}
	return x;
}






/**
 *测试矩阵运算成功
 */
void test_Mat(){
    vector<vector<Double> > data,tdata,res,Ax;
    //data = getdata();
	Double x[] = {2,1,-1,2,1,0,1,-1,1};
	
    for(int i=0;i<3;i++){
        vector<Double> tmp;
        data.push_back(tmp);
        for(int j=0;j<3;j++){
            data[i].push_back(x[i*3+j]);
        }
    }
	
	/**
    tdata = t(data);

    for(int i=0;i<tdata.size();i++){
        for(int j=0;j<tdata[0].size();j++){
            cout<<tdata[i][j]<<" ";
        }
        cout<<endl;
    }

    res = mulMat(tdata,data);
    for(int i=0;i<res.size();i++){
        for(int j=0;j<res[0].size();j++){
            cout<<res[i][j]<<" ";
        }
        cout<<endl;
    }
	
	cout<<det(data)<<endl;
	*/
	
	//Ax = inv(data);
	cout<<det(data)<<endl;
	cout<<"逆矩阵:"<<endl;
	for(int i=0;i<Ax.size();i++){
        for(int j=0;j<Ax[0].size();j++){
            cout<<Ax[i][j]<<" ";
        }
        cout<<endl;
    }
}

/**
int main()
{

    vector<Double> data;


    test_Mat();
    return 0;
}

**/



