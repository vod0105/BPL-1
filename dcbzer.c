#include <stdio.h>
#include <math.h>

// B la nghich dao cua ma tran he so cua he phuong trinh hai an (alpha1*t1) va (alpha2*t2)
// D la ma tran cot cua he so tu do cua he phuong trinh hai an (alpha1*t1) va (alpha2*t2)
double B[2][2], D[2][3];

// u tuong ung voi diem a[i] nam tren duong cong
double t(int n, int i);

// tinh B[i,n](u)_ da thuc Bernstein
double B0(double u);
double B1(double u);
double B2(double u);
double B3(double u);

double distanse(double a[100][3], int i);         // khoang cach tu diem den tam 0
void sort(double a[100][3], int n);               // sap xep cac diem theo thu tu
void Swap(double a[100][3], int n, int i, int j); // doi dong i cho j

// tinh tong cac B[i,n] cua tat ca cac u
double tongB0(double a[100][3], int n);
double tongB1(double a[100][3], int n);
double tongB2(double a[100][3], int n);
double tongB3(double a[100][3], int n);

double tongB1B1(double a[100][3], int n);
double tongB1B2(double a[100][3], int n);
double tongB2B2(double a[100][3], int n);

// tinh ma tran B va ma tran D
void MTB(double a[100][3], int n);
void MTD(double a[100][3], int n);

void output(double a[100][3], int n);
void readfile(double a[100][3], int *n);
void writefile(double P[4][3]);

int main()
{
    // m la lua chon
    // n la so cac diem nam tren duong cong
    // a la toa do cac diem nam tren duong cong
    int m, n;
    double a[100][3];
    printf("_____________________________________________________________________ \n");
    printf("|                     NOI SUY DUONG CONG BEZIER                      |\n");
    printf("|     Nguoi HD: PGS. TS. NGUYEN TAN KHOI                             |\n");
    printf("|     Sinh vien TH:                                                  |\n");
    printf("|     VO VAN DUC- 21T_DT2                                            |\n");
    printf("|     NGUYEN THI THU THAO- 21T_DT2                                   |\n");
    printf("|____________________________________________________________________| \n\n\n\n");
    do
    {
        printf("_________----MENU----___________\n");
        printf("AN 1:   Lay  du  lieu  tu   file\n");
        printf("AN 2:   Nhap du lieu tu ban phim\n");
        printf("AN 3:   Thoat\n\n");
        printf("Nhap lua chon cua ban: ");
        scanf("%d", &m);

        // nhap du lieu......................
        switch (m)
        {
        // lay du lieu tu file txt...............
        case 1:
            readfile(a,&n);
            break;

        // nhap du lieu truc tiep.....................
        case 2:
            printf("Nhap so diem tren duong cong: ");
            scanf("%d", &n);

            printf("nhap toa do cac diem tren duong cong:\n");
            for (int i = 0; i < n; i++)
            {
                printf("Nhap toa do diem thu %d:\n", i);
                for (int j = 0; j < 3; j++)
                    scanf("%lf", &a[i][j]);
            }

            printf("Cac diem tren duong cong la:\n");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 3; j++)
                    printf("%.2lf\t", a[i][j]);
                printf("\n");
            }
            break;

        case 3:
            return 0;
        }

        // xuat ket qua.................................
        sort(a, n);
        output(a, n);

    } while (1);

    return 0;
}

double t(int n, int i)
{
    return (double)i / (n - 1);
}

//....................................
double B0(double u)
{
    return (1 - u) * (1 - u) * (1 - u);
}

double B1(double u)
{
    return 3 * (1 - u) * (1 - u) * u;
}

double B2(double u)
{
    return 3 * (1 - u) * u * u;
}

double B3(double u)
{
    return u * u * u;
}

// doi dong............
void Swap(double a[100][3], int n, int i, int j)
{
    double tmp;
    for (int e = 0; e < 3; e++)
    {
        tmp = a[i][e];
        a[i][e] = a[j][e];
        a[j][e] = tmp;
    }
}

// khoang cach cua diem den tam 0...
double distance(double a[100][3], int i)
{
    double e = sqrt(a[i][0] * a[i][0] + a[i][1] * a[i][1] + a[i][2] * a[i][2]);
    return e;
}

void sort(double a[100][3], int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {

            double e = distance(a, i);
            double f = distance(a, j);
            if (e > f)
                Swap(a, n, i, j);
        }
    }
}

//.........................
double tongB0(double a[100][3], int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        double u = t(n, i);
        sum += B0(u);
    }
    return sum;
}

double tongB1(double a[100][3], int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        double u = t(n, i);
        sum += B1(u);
    }
    return sum;
}

double tongB2(double a[100][3], int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        double u = t(n, i);
        sum += B2(u);
    }
    return sum;
}

double tongB3(double a[100][3], int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        double u = t(n, i);
        sum += B3(u);
    }
    return sum;
}

//.........................
double tongB1B1(double a[100][3], int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        double u = t(n, i);
        sum += B1(u) * B1(u);
    }
    return sum;
}

double tongB2B2(double a[100][3], int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        double u = t(n, i);
        sum += B2(u) * B2(u);
    }
    return sum;
}

double tongB1B2(double a[100][3], int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        double u = t(n, i);
        sum += B1(u) * B2(u);
    }
    return sum;
}

// Ma tran he so cua he phuong trinh ...........................
void MTB(double a[100][3], int n)
{
    double sum1 = tongB1B1(a, n);
    double sum2 = tongB2B2(a, n);
    double sum3 = tongB1B2(a, n);
    double sum = sum1 * sum2 - sum3 * sum3;
    B[0][0] = sum2 / sum;
    B[1][0] = -sum3 / sum;
    B[0][1] = -sum3 / sum;
    B[1][1] = sum1 / sum;
}

// Ma tran ket qua cua he phuong trinh............................
void MTD(double a[100][3], int n)
{
    for (int k = 0; k < 3; k++)
    {
        double sum1 = 0;
        for (int i = 0; i < n; i++)
            sum1 += (a[i][k] - (a[0][k] * B0(t(n, i)) + a[0][k] * B1(t(n, i)) + a[n - 1][k] * B2(t(n, i)) + a[n - 1][k] * B3(t(n, i)))) * B1(t(n, i));
        D[0][k] = sum1;

        double sum2 = 0;
        for (int i = 0; i < n; i++)
            sum2 += (a[i][k] - (a[0][k] * B0(t(n, i)) + a[0][k] * B1(t(n, i)) + a[n - 1][k] * B2(t(n, i)) + a[n - 1][k] * B3(t(n, i)))) * B2(t(n, i));
        D[1][k] = sum2;
    }
}

//...............................
void output(double a[100][3], int n)
{
    MTB(a, n);
    MTD(a, n);

    // ma tran cot cua an (alpha*t)..................................
    double alpha_t[2][3];
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++)
            alpha_t[i][j] = B[i][0] * D[0][j] + B[i][1] * D[1][j];

    // ma tran cua cac diem dieu khien..............................
    double P[4][3];
    for (int i = 0; i < 3; i++)
        P[0][i] = a[0][i];
    for (int i = 0; i < 3; i++)
        P[3][i] = a[n - 1][i];
    for (int i = 0; i < 3; i++)
        P[1][i] = alpha_t[0][i] + P[0][i];
    for (int i = 0; i < 3; i++)
        P[2][i] = alpha_t[1][i] + P[3][i];

    // xuat cac diem dieu khien.........................................
    printf("Cac diem dieu khien la : \n");
    for (int i = 0; i < 4; i++)
        printf("P[%d]= (%10.2lf ,%10.2lf ,%10.2lf)\n", i, P[i][0], P[i][1], P[i][2]);

    //ghi ket qua vao file
    int tmp;
    printf("Ban co muon ghi ket qua vao file ko :\n");
    printf("1: YES\n");
    printf("2: NO\n");
    printf("nhap lua chon:\t");
    scanf("%d",&tmp);
    if(tmp==1) writefile(P);
     printf ("Nhiem vu hoan thanh\n");           
}

//doc file

void readfile(double a[100][3],int *n){
           char file[100];
            FILE *f;
            printf("Nhap ten file :");fflush(stdin);
            gets(file);
            f = fopen(file, "r");
            fscanf(f, "%d ", n);
            for (int i = 0; i < *n; i++)
                for (int j = 0; j < 3; j++)
                    fscanf(f, "%lf ", &a[i][j]);
            fclose(f);        
            printf("Cac diem tren duong cong la:\n");
            for (int i = 0; i < *n; i++)
            {
                for (int j = 0; j < 3; j++)
                    printf("%.2lf\t", a[i][j]);
                printf("\n");
            }
}

void writefile(double P[4][3]){
           char file[100];
           FILE *f;
           printf("Nhap file can ghi ket qua vao: ");
           fflush(stdin);
           gets(file);
           f=fopen(file,"w");
           fprintf(f,"\n cac diem dieu khien la:\n ");
           for (int i=0;i<4;i++){
           fprintf(f,"P[%d]= (%10.2lf ,%10.2lf ,%10.2lf)\n", i, P[i][0], P[i][1], P[i][2]);
           }
           fclose(f);
}