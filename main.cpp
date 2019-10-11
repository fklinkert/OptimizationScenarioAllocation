#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <gurobi_c++.h>
using namespace std;

const int ROW = 756;
const int COL = 487;

void getYX(double** price, string line, char delimeter, int row){
    // split the whole line by delimeter and store them into array of doubles
    int i = 0;
    int start = 0;
    int count = 0;
    
    while(i < line.length()){
        if(line[i] == delimeter){
            string str = line.substr(start, i-start);
            start = i+1;
            double val = 0;
            if(str != "") val = stod(str);
            if(count == 0) price[row][0] = val;
            if(count > 1 &&  row < ROW) {
                if(row == 30){
                }
                price[row][count-1]= val;
            }
            count ++;
        }
        i ++;
    }
    if(row < ROW){
        string str = "";
        if(line.length()-1 > start)
            str = line.substr(start, line.length()-start);
        double val = 0;
        if(str != "") val = stod(str);
        price[row][count-1]= val;
    }
}

void getData(double** price, string filename){
    // store data into index and price
    ifstream file(filename);
    
    string line = "";
    
    int i = 0;
    while (getline(file, line)){
        if(i > 0) {
            getYX(price, line, ',', i-1);
        }
        i ++;
        
    }
    file.close();
}

void cov_with_risk_free_asset(double** price, double** cov, int p){
    
    double means[COL];
    
    for(int j = 0; j < COL; j++){
        means[j] = 0;
        for(int i = int(ROW/4)*p; i < int(ROW/4)*p+int(ROW/4); i ++){
            means[j] += price[i][j];
        }
        means[j] = 4*means[j]/ROW;
    }
    
    double sum;
    
    for(int j = 0; j< COL; j++){
        for(int k = 0; k < COL; k++){
            sum = 0;
            for(int i =int(ROW/4)*p ; i < int(ROW/4)+int(ROW/4)*p; i++)
                sum += (price[i][j] - means[j]) * (price[i][k] - means[k]);
            cov[j][k] = sum/(int(ROW/4));
        }
    }
    
}


void matrixMult(double** matrix1, double** matrix2, double** result, int rowFirst, int columnSecond, int columnFirst){
    
    for(int i = 0; i < rowFirst; ++i)
    {
        for(int j = 0; j < columnSecond; ++j)
        {
            result[i][j] = 0;
        }
    }
    
    // Multiplying matrix1 and matrix2 and storing in array mult.
    for(int i = 0; i < rowFirst; ++i)
    {
        for(int j = 0; j < columnSecond; ++j)
        {
            for(int k=0; k<columnFirst; ++k)
            {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    
}

void MovingAverage(double** data, double** MovingAvg){
    double alpha = 0.5;
    
    for(int i = 0; i < ROW; i++){
        for(int j = 0; j < COL; j++){
            if(i == 0)
                MovingAvg[i][j] = data[i][j];
            else
                MovingAvg[i][j] = data[i][j]*alpha + (1-alpha)*MovingAvg[i-1][j];
        }
    }
}

int main(int argc, const char * argv[])
{
    
     // data reading
     if (argc != 2){
     cerr << "usage: ./a.out <filename>\n";
     exit(1);
     }
     
    cout << fixed;
    cout << setprecision(8);
    
    // read file and store
    string filename = argv[1];
    
    double** data = new double*[ROW];
    for(int i = 0; i < ROW; ++i)
        data[i] = new double[COL];
    
    getData(data, filename);
    
    
    double* TenDaysPerformance = new double[COL];
    double performance = 0;
    
    for(int j = 0; j < COL; j++){
        if(data[0][j] != 0)
            TenDaysPerformance[j] = (data[9][j] - data[0][j])/data[0][j];
        else
            TenDaysPerformance[j] = 0;
    }
    
    
    
    double** covariance = new double*[COL];
    for (int i = 0; i < COL; i++)
        covariance[i] = new double[COL];
    
    double** moving_avg = new double*[ROW];
    for(int i = 0; i < ROW; i++){
        moving_avg[i] = new double[COL];
    }
    
    MovingAverage(data, moving_avg);
    
    double** p = new double*[COL];
    for(int i =0; i < COL; i++)
        p[i] = new double[COL];
    
    for(int i = 0; i < COL; i++){
        for(int j = 0; j< COL; j++){
            if(i==j)
                p[i][j] = 1/data[0][i];
            else
                p[i][j] = 0;
        }
    }
    
    
    
    int retcode = 0;
    double B = 1000;
    
    char** mynames = new char*[COL+2];
    for(int j = 0; j < COL+1; j++){
        mynames[j] = new char[4];
        sprintf(mynames[j],"z%d", j);
    }
    mynames[COL+1] = new char[4];
    sprintf(mynames[COL+1],"std");
    
    double* obj = new double[COL+1];
    double* ub = new double[COL+2];
    double* lb = new double[COL+2];
    double* z = new double[COL+2];
    
    for(int i = 0; i < COL+1; i++){
        obj[i] = 1.0;
        ub[i] = data[0][i]/B;
        lb[i] = (-data[0][i]/B);
    }
    
    ub[COL+1] = 10000;
    lb[COL+1] = 0;
    
    /** we need to have a couple of auxiliary arrays **/
    int* cind = new int[COL+2];
    double* cval = new double[COL+2];
    
    int numnonz;
    double rhs;
    char sense;
    
    double sig = 0;
    double vega = 0;
    
    
    double** MovingAvg_period = new double*[int(ROW/4)];
    for(int i =0; i < int(ROW/4); i++)
        MovingAvg_period[i] = new double[COL];
    
    double** mult = new double*[COL];
    double** cov_z = new double*[COL];
    for(int i =0; i < COL; i++){
        mult[i] = new double[COL];
        cov_z[i] = new double[COL];
        
    }
    
    double** covariance_period = new double*[COL];
    for(int i =0; i < COL; i++)
        covariance_period[i] = new double[COL];
    
    int quadratic = COL*COL+1;
    
    int* qrow = new int[quadratic];
    int* qcol = new int[quadratic];
    double* qval = new double[quadratic];
    
    for(int i = 0; i < COL+1; i++){
        for(int j = 0; j < COL+1; j++){
            if(i!=0 && j!=0){
                qrow[COL*(i-1)+(j-1)] = j;
                qcol[COL*(i-1)+(j-1)] = i;
            }
        }
    }
    
    qrow[quadratic-1] = COL+1;
    qcol[quadratic-1] = COL+1;
    qval[quadratic-1] = -1;
    
    
    
    for(int i = 0; i < 6; i++){
        for(int j = 1; j < 4; j++){
            GRBenv   *env = NULL;
            GRBmodel *model = NULL;
            
            retcode = GRBloadenv(&env, "arbitrage.log");
            if (retcode) return retcode;
            
            retcode = GRBnewmodel(env, &model, "arbitrage", COL+2, NULL, NULL, NULL, NULL, NULL);
            if (retcode) return retcode;
            
            /* initialize variables */
            for(int j = 0; j < COL+2; j++){
                retcode = GRBsetstrattrelement(model, "VarName", j, mynames[j]);
                if (retcode) return retcode;
                
                retcode = GRBsetdblattrelement(model, "Obj", j, obj[j]);
                if (retcode) return retcode;
                
                retcode = GRBsetdblattrelement(model, "LB", j, lb[j]);
                if (retcode) return retcode;
                
                retcode = GRBsetdblattrelement(model, "UB", j, ub[j]);
                if (retcode) return retcode;
            }
            
            sig = 0.2 * i;
            vega = 0.5 * j;
            
            /** first constraint is: budget constraint **/
            
            for(int i = 0; i < COL+2; i++){
                cind[i] = i;
                cval[i] = 1;
            }
            cval[COL+1] = -vega;
            
            numnonz = COL+1;
            rhs = (1 - 2 * sig);
            sense = GRB_EQUAL;
            
            retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "budget_constraint");
            if (retcode) return retcode;
            
            /** second constraint is: chance constraint period T/4**/
            for(int i = 0; i < int(ROW/4); i++)
                for(int j = 0; j < COL; j++)
                    MovingAvg_period[i][j] = moving_avg[i][j];
            
            for(int i = 0; i < COL+1; i++){
                if(i ==0)
                    cval[i] = pow(1.0+0.00013,189);
                else
                    cval[i] = MovingAvg_period[188][i-1]/MovingAvg_period[0][i-1];
            }
            
            numnonz = COL+2;
            rhs = 0.0;
            sense = GRB_GREATER_EQUAL;
            
            retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "chance_constraint_1");
            if (retcode) return retcode;
            
            
            /** third constraint is: variance constraint period T/4 **/
            
            cov_with_risk_free_asset(MovingAvg_period, covariance_period, 0);
            
            matrixMult(p, covariance_period, mult, COL, COL, COL);
            matrixMult(mult, p, cov_z, COL, COL, COL);
            
            
            for(int i = 0; i < COL; i++){
                for(int j = 0; j < COL; j++){
                    qval[COL*i+j] = cov_z[i][j];
                }
            }
            
            retcode = GRBaddqconstr(model, 0, NULL, NULL, quadratic, qrow, qcol, qval,GRB_LESS_EQUAL, rhs, "variance_constraint_1");
            if (retcode) return retcode;
            
    
            /** second constraint is: chance constraint period T/2**/
            for(int i = 0; i < int(ROW/4); i++)
                for(int j = 0; j < COL; j++)
                    MovingAvg_period[i][j] = moving_avg[i+int(ROW/4)][j];
            
            for(int i = 0; i < COL+1; i++){
                if(i ==0)
                    cval[i] = pow(1.0+0.00013,2*189);
                else
                    cval[i] = MovingAvg_period[188][i-1]/MovingAvg_period[0][i-1];
            }
            
            numnonz = COL+2;
            rhs = 0.0;
            sense = GRB_GREATER_EQUAL;
            
            retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "chance_constraint_2");
            if (retcode) return retcode;
            
            
            /** third constraint is: variance constraint period T/2 **/
            
            cov_with_risk_free_asset(MovingAvg_period, covariance_period, 0);
            
            
            matrixMult(p, covariance_period, mult, COL, COL, COL);
            matrixMult(mult, p, cov_z, COL, COL, COL);
            
            
            for(int i = 0; i < COL; i++){
                for(int j = 0; j < COL; j++){
                    qval[COL*i+j] = cov_z[i][j];
                }
            }
            
            retcode = GRBaddqconstr(model, 0, NULL, NULL, quadratic, qrow, qcol, qval,GRB_LESS_EQUAL, rhs, "variance_constraint_2");
            if (retcode) return retcode;
            
            /** second constraint is: chance constraint period 3T/4**/
            for(int i = 0; i < int(ROW/4); i++)
                for(int j = 0; j < COL; j++)
                    MovingAvg_period[i][j] = moving_avg[i+2*int(ROW/4)][j];
            
            for(int i = 0; i < COL+1; i++){
                if(i ==0)
                    cval[i] = pow(1.0+0.00013,3*189);
                else
                    cval[i] = MovingAvg_period[188][i-1]/MovingAvg_period[0][i-1];
            }
            
            numnonz = COL+2;
            rhs = 0.0;
            sense = GRB_GREATER_EQUAL;
            
            retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "chance_constraint_3");
            if (retcode) return retcode;
            
            
            /** third constraint is: variance constraint period 3T/4 **/
            
            cov_with_risk_free_asset(MovingAvg_period, covariance_period, 0);
            
            
            matrixMult(p, covariance_period, mult, COL, COL, COL);
            matrixMult(mult, p, cov_z, COL, COL, COL);
            
            
            for(int i = 0; i < COL; i++){
                for(int j = 0; j < COL; j++){
                    qval[COL*i+j] = cov_z[i][j];
                }
            }
            

            
            retcode = GRBaddqconstr(model, 0, NULL, NULL, quadratic, qrow, qcol, qval,GRB_LESS_EQUAL, rhs, "variance_constraint_3");
            if (retcode) return retcode;
            
            
            
            retcode = GRBupdatemodel(model);
            if (retcode) return retcode;
            
            /** optional: write the problem **/
            retcode = GRBwrite(model, "arbitrage.lp");
            if (retcode) return retcode;
            
            retcode = GRBoptimize(model);
            if (retcode) return retcode;
            
            /** get solution **/
            retcode = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, COL+2, z);
            
            /** now let's see the values **/
            for(int j = 0; j < COL+2; j++){
                printf("%s = %g\n", mynames[j], z[j]);
            }
            
            performance = 0;
            
            for(int j = 0;j< COL ;j++)
                performance+= TenDaysPerformance[j]*z[j+1];
            performance+= z[0]*pow(1.0+0.00013,10);
            cout << "\nThe performance of the portfolio after ten days equals " << performance <<
                "% of the budget\n";
            
            GRBfreeenv(env);
        }
    }
    return retcode;
}

