#include <iostream>
#include <fstream>
#include<string>
#include<math.h>
#include<stdlib.h>
using namespace std;
const double g = 1.4;  // gamma
typedef double (*Fun) (double);
double Nlsolve(Fun f, Fun df, double itl){
    double ret = itl;
    double err = 1.;
    int n = 0;
    while(err>1e-12){
        err = f(ret) / df(ret);
        ret -= err;
        err = abs(err);
        n++;
        cout << "Iteration = " << n << " ; err = " << err << " ;" << endl;
    }
    return ret;
}
// 左侧状态
double p1 ;
double rho1 ;
double u1 ;
double c1;

// 右侧状态
double p2 ;
double rho2 ;
double u2 ;
double c2;

// 中间状态
double ustar;
double pstar;
double rhoL;
double cL;
double rhoR;
double cR;

bool RS;  //用于判定右侧激波
double Z1;
bool LS;  //用于判断左侧激波
double Z2;

bool isVac;  //用于判断中间是否真空
void SetLetfState(double p,double rho,double u){
    p1 = p;
    rho1 = rho;
    u1 = u;
}
void SetRightfState(double p,double rho,double u){
    p2 = p;
    rho2 = rho;
    u2 = u;
}
void SetState(double _p1,double _rho1,double _u1,double _p2,double _rho2,double _u2){
    p1 = _p1;
    rho1 = _rho1;
    u1 = _u1;
    p2 = _p2;
    rho2 = _rho2;
    u2 = _u2;
}


double f(double p,double pi,double rhoi){
    double ret,ci;
    ci = sqrt(g * pi/rhoi);
    if(p>pi){
        ret = (p-pi)/(rhoi*ci*sqrt(((g+1)*p/pi+(g-1))*0.5/g));
    }
    else{
        ret = 2*ci/(g-1)*(pow(p/pi,(g-1)*0.5/g)-1);
    }
    return ret;
}
double df(double p,double pi,double rhoi){
    double ret,ci;
    ci = sqrt(g * pi/rhoi);
    if(p>pi){
        ret = ((g+1)*p/pi + (3*g-1))/(4*g*rhoi*ci*pow(((g+1)*p/pi+(g-1))/(2*g),1.5));
    }
    else{
        ret = ci/(g*pi)*pow(p/pi,-(g+1)/(2*g));
    }
    return ret;
}
double F(double p){
    return (f(p,p1,rho1) + f(p,p2,rho2) - (u1-u2));
}
double dF(double p){
    return (df(p,p1,rho1) + df(p,p2,rho2));
}

//double f_test(double x){
//    return (sin(x) + x - 0.5);
//}
//double df_test(double x){
//    return (cos(x) + 1);
//}

double rho_rs(double pr){ //激波前后的密度比 rho ratio of shock
    double k = (g+1)/(g-1);
    return (k*pr+1)/(k+pr);
}
double rho_re(double pr){
    return pow(pr,1/g);
}
void deteparams(){
    c1 = sqrt(g*p1/rho1);
    c2 = sqrt(g*p2/rho2);

    cout << F(p1) <<endl;
    if(F(0.)>0.){ //中间真空
        isVac = true;
        pstar = 0.;
        rhoL = 0.;
        rhoR = 0.;
        ustar = 0.;
        cL = 0.;
        cR = 0.;
    }
    else {
        isVac = false;
        // 计算p*
        
        pstar = Nlsolve(F,dF,1e-8); //计算中间的压强
        ustar = 0.5*(u1+u2 + f(pstar,p2,rho2)- f(pstar,p1,rho1));  //计算中间的速度
        cout << pstar <<endl;
        //cout << ustar <<endl;
        
        LS = pstar>p1;
        RS = pstar>p2;
        //计算右侧
        if(RS){
            rhoR = rho2*rho_rs(pstar/p2); // 右侧激波
            Z2 = (rhoR*ustar - rho2*u2)/(rhoR-rho2);
        }
        else{
            rhoR = rho2*rho_re(pstar/p2); //右侧膨胀波
            Z2 = nan("");
        }
        cR = sqrt(g*pstar/rhoR);

        //计算左侧
        if(LS){
            rhoL = rho1*rho_rs(pstar/p1);  //左侧激波
            Z1 = (rhoL*ustar - rho1*u1)/(rhoL-rho1);
        }
        else{
            rhoL = rho1*rho_re(pstar/p1); //左侧膨胀波
            Z1 = nan("");
        }
        cL = sqrt(g*pstar/rhoL);
    }
}
void getState(double t,double *x,double *p,double *rho,double *u){ //计算t时刻，x位置的状态
    if(isVac){ //中间真空
        if(*x<(u1+2*c1/(g-1))*t){
            if(*x<(u1-c1)*t){
                *u = u1;*p=p1;*rho = rho1;
            }
            else{
                double c = ((g-1)*(u1 - *x/t)+2*c1)/(g+1);
                *u = c + *x/t;
                *p = p1*pow(c/c1,2*g/(g-1));
                *rho = g*(*p)/(c*c);
            }
        }
        else if (*x>(u2-2*c2/(g-1))*t){
            if(*x>(u2+c2)*t){
                *u = u2;*p=p2;*rho = rho2;
            }
            else{
                double c = ((g-1)*(*x/t-u2)+2*c2)/(g+1);
                *u = *x/t - c;
                *p = p2*pow(c/c2,2*g/(g-1));
                *rho = g*(*p)/(c*c);
            }
        }
        else{
            *u = ustar;
            *p = pstar;
            *rho = 0.0;
        }
    }
    else{
        if(*x>ustar*t){ //计算右侧
            if(RS){  //右侧激波
                if(*x<Z2*t){
                    *u = ustar;*p=pstar;*rho = rhoR;
                }else{
                    *u = u2;*p=p2;*rho = rho2;
                }
            }
            else{  //右侧膨胀波
                if(*x<(ustar+cR)*t){
                    *u = ustar;*p=pstar;*rho = rhoR;
                }
                else if(*x<(u2+c2)*t){
                    double c = ((g-1)*(*x/t-u2)+2*c2)/(g+1);
                    *u = *x/t - c;
                    *p = p2*pow(c/c2,2*g/(g-1));
                    *rho = g*(*p)/(c*c);
                }
                else{
                    *u = u2;*p=p2;*rho = rho2;
                }
            }
        }       
        else{ //计算左侧
            if(LS){ //左侧激波
                if(*x>Z1*t){
                    *u = ustar;*p=pstar;*rho = rhoL;
                }else{
                    *u = u1;*p=p1;*rho = rho1;
                }   
            }
            else{ //左侧膨胀波
                if(*x<(u1-c1)*t){
                    *p = p1;*rho = rho1;*u = u1;
                }else if(*x<(ustar-cL)*t){
                    double c = ((g-1)*(u1 - *x/t)+2*c1)/(g+1);
                    *u = c + *x/t;
                    *p = p1*pow(c/c1,2*g/(g-1));
                    *rho = g*(*p)/(c*c);
                }else {
                    *u = ustar;*p=pstar;*rho = rhoL;
                }
            }

        }
    }
}
void getdata(double t,double *x,double *p,double *rho,double *u,int size){
    for(int i =0;i<size;i++){
        getState(t,x,p,rho,u);
        x++;p++;rho++;u++;
    }
}
void outdata_at_t(string filename,double *x,double *p,double *rho,double *u,int size){
    //把t时刻数据输出到文件
    ofstream csvfile;
    csvfile.open(filename, ios::out | ios::trunc);
    csvfile << "x" << "," << "p" << "," << "rho" << "," << "u" << endl;
    for(int i=0;i<size;i++){
        csvfile << *x << "," << *p << "," << *rho << "," << *u << endl;
        x++;p++;rho++;u++;
    }  
};
void pArr2f(string filename, double *p,int tsize,int xsize){ //把二维数组输出到文件
    ofstream csvfile;
    csvfile.open(filename,ios::out | ios::trunc);
    for(int j=0;j<xsize;j++){
        for(int i=0;i<tsize;i++){
            csvfile << p[j + xsize*i] <<",";
        }
        csvfile << endl;
    }
}
void outdata_span_t(double *t,double *x,double *p,double *rho,double *u,int tsize,int xsize){
    //输出时间段的数据
    pArr2f("t.csv",t,1,tsize);
    pArr2f("x.csv",x,1,xsize);
    pArr2f("p.csv",p,tsize,xsize);
    pArr2f("rho.csv",rho,tsize,xsize);
    pArr2f("u.csv",u,tsize,xsize);
}
int main(){

    SetState(1,1,-0.,0.1,0.125,0.);
    deteparams();
    if(LS){
        cout<<"左行激波"<<endl;
    }else{
        cout<<"左行膨胀波"<<endl;
    }
    if(RS){
        cout<<"右行激波"<<endl;
    }else{
        cout<<"右行膨胀波"<<endl;
    }
    if(isVac) cout << "中间真空" <<endl;
    cout << ustar<<endl;
    cout << Z1 <<endl;
    cout <<  Z2 <<endl;
    double t = 1.4;
    const int n = 401;
    const double L = 3*ceil(t*max(u2+c2,max(abs(u1-c1),max(abs(Z1),abs(Z2)))));
    cout << L <<endl;
    double x[n],p[n],rho[n],u[n];

    double d = L/(n-1);
    for(int i=0;i<n;i++){
        x[i] = d*i - L*0.5;
    }
    getdata(t,x,p,rho,u,n);  //获取t时刻的数据
    outdata_at_t("test_t=1.4.csv",x,p,rho,u,n);  //输出t=1.4s的数据

    double dt = 0.05;
    const int nt = 41;
    double t_end = (nt-1)*dt;
    const double Ls = 3*ceil(t_end*max(u2+c2,max(abs(u1-c1),max(abs(Z1),abs(Z2)))));
    const int nx = 401;
    double dx = Ls/(nx-1);

    double ts[nt],xs[nx];
    double ps[nt][nx],rhos[nt][nx],us[nt][nx];
    for(int i=0;i<nx;i++){
        xs[i] = dx*i - Ls*0.5;
    }
    for(int i=0;i<nt;i++){
        ts[i] = i*dt;
    }
    
    for(int i=0;i<nt;i++){
        getdata(ts[i],xs,ps[i],rhos[i],us[i],nx);
    }
    outdata_span_t(ts,xs,ps[0],rhos[0],us[0],nt,nx); //输出各时刻的数据

    
    // 测试用
    //cout << f(1.0) << endl;
    //cout << df(1.0) << endl;
    //double sol = Nlsolve(f,df,0.5);
    //cout << sol << endl;
    //cout << f(sol);
    system("pause");

    return 0;
}