
#pragma once


class Interpolation //�ߴ��������ղ�ֵ
{
public:
    static double lgr(double const *lgrx, double const *lgry, int lgrn, double lgrt);
    static double lg3(double const *x, double const *y, int n, double t);
    static double Linear_interpolation(double const *x, double const*y, int n, double t);
    //   double Linear_interpolation1(double *x, double *y, int n, double t);
};


