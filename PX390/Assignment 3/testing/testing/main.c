//
//  main.c
//  testing
//
//  Created by Ajan Sittampalam on 20/10/2020.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double diff(double x1, double x2, double y1, double y2, double z1, double z2){
    double d;
    d = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    d = sqrt(d);
    return d;
}

void en_pot(double *posx, double *posy, double *posz,  long ncharges, double *res){
    double p = 0;
    for (int i=0;i<ncharges;i++){
        for (int j=0;j<ncharges;j++){
            if (i==j){
                p = p +0;
            }
            else{
                p = p+ 1/diff(posx[i],posx[j],posy[i],posy[j],posz[i],posz[j]);
            }
        }
    }
    *res = 0.5*p;
    return;
}
int main(void) {
    long charges =2;
    double x[]={1,2};
    double y[]={5,6};
    double z[]={3,4};
    double res;
    en_pot(x, y, z, charges , &res);
    printf("%lf ",res);
}
