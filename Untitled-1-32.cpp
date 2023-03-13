#include <iostream>
#include <fstream>
using namespace std;
#include<string>
#include<math.h>
#include<stdlib.h>
typedef float (*Fun) (float);
float g = 1.4;
float p1  = 1.;
float rho1 = 1.;
float p2 = 0.1;
float rho2 = 0.125;
float u1 = 0.0;
float u2 = 0.0;
float ustar;
float pstar;

float f(float p,float pi,float rhoi){
    float ret,ci;
    ci = sqrt(g * pi/rhoi);
    if(p>pi){
        ret = (p-pi)/(rhoi*ci*sqrt(((g+1)*p/pi+(g-1))*0.5/g));
    }
    else{
        ret = 2*ci/(g-1)*(pow(p/pi,(g-1)*0.5/g)-1);
    }
    return ret;
}
float df(float p,float pi,float rhoi){
    float ret,ci;
    ci = sqrt(g * pi/rhoi);
    if(p>pi){
        ret = ((g+1)*p/pi + (3*g-1))/(4*g*rhoi*ci*pow((g+1)*p/pi+(g-1),1.5));
    }
    else{
        ret = ci/(g*pi)*pow(p/pi,-(g+1)/(2*g));
    }
    return ret;
}
float F(float p){
    return (f(p,p1,rho1) + f(p,p2,rho2) - (u1-u2));
}
float dF(float p){
    return (df(p,p1,rho1) + df(p,p2,rho2));
}
//float f_test(float x){
//    return (sin(x) + x - 0.5);
//}
//float df_test(float x){
//    return (cos(x) + 1);
//}
float Nlsolve(Fun f, Fun df, float itl){
    float ret = itl;
    float err = 1.;
    int n = 0;
    while(err>1e-6){
        err = f(ret) / df(ret);
        ret -= err;
        err = abs(err);
        n++;
        cout << "Iteration = " << n << " ; err = " << err << " ;" << endl;
    }
    return ret;
}

int main(){
    // 测试用
    //cout << f(1.0) << endl;
    //cout << df(1.0) << endl;
    //float sol = Nlsolve(f,df,0.5);
    //cout << sol << endl;
    //cout << f(sol);

    // 计算p⭐
    pstar = Nlsolve(F,dF,0.5); //计算中间的压强
    ustar = u1 - f(pstar,p1,rho1);  //计算中间的速度
    //cout << pstar <<endl;
    //cout << F(pstar)<<endl;
    //cout << ustar <<endl;
    ofstream csvfile;
    csvfile.open("test.csv", ios::out | ios::trunc);
    for(int i=1;i<=200;i++){
        csvfile << 0.05*i << "," << F(0.05*i) << endl;
    }  
    return 0;
}