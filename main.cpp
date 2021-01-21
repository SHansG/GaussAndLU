#include<iostream>
#include<iomanip>
#include<math.h>

using namespace std;

//Utility


float** get_matrix(int size)
{
    float** matrix = new float* [size];

    for (int i = 0; i < size; i++)
    {
        matrix[i] = new float[size];
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[i][j] = 0.0f;
        }
    }
    return matrix;
}

float* get_tab(int size)
{
    float* tab = new float[size];
    for (int i = 0; i < size; i++)
    {
        tab[i] = 0.0;
    }
    return tab;
}

void show_matrix(float** matrix, int& size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
            printf("%6.1f ", matrix[i][j]);
        printf("\n");
    }
}

void show_coefficient_matrix(float** matrix, int& size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j <= size; j++)
            printf("%6.1f ", matrix[i][j]);
        printf("\n");
    }
}

//Gauss

void forward_elim(float** matrix, int& size)
{
    float tmp;
    for (int k = 0; k < size - 1; k++)
    {
        for (int i = k; i < size - 1; i++)
        {
            tmp = (matrix[i + 1][k] / matrix[k][k]);

            for (int j = 0; j <= size; j++)
            {
                matrix[i + 1][j] -= tmp * matrix[k][j];
            }
        }
    }
}

void backward_substitution(float** matrix, float* result_tab, int& size)
{
    float tmp;
    for (int i = size - 1; i >= 0; i--)
    {
        tmp = 0;
        for (int j = i; j <= size - 1; j++)
        {
            tmp = tmp + matrix[i][j] * result_tab[j];
        }
        result_tab[i] = (matrix[i][size] - tmp) / matrix[i][i];
    }
}

// LU

void upper_triangular_matrix(float** A, float** upper, float** lower, int& size ,int& row)
{
    float tmp;
    for (int i = row; i < size; i++)
    {
        tmp = 0;
        for (int j = 0; j < row; j++)
        {
            tmp += upper[j][i] * lower[row][j];
        }
        upper[row][i] = A[row][i] - tmp;
    }
}

void lower_triangular_matrix(float** A, float** upper, float** lower, int& size, int& col)
{
    float tmp;
    for (int i = col + 1; i < size; i++)
    {
        tmp = 0;
        for (int j = 0; j <= i - 1; j++)
        {
            tmp += upper[j][col] * lower[i][j];
        }
        lower[i][col] = (A[i][col] - tmp) / upper[col][col];
    }
}

//LY = B
void LY_B(float** lower, float* B, float* Y, int& size)
{
    float tmp;
    for (int i = 0; i < size; i++)
    {
        tmp = 0;
        for (int j = 0; j <= i - 1; j++)
        {
            tmp += lower[i][j] * Y[j];
        }
        Y[i] = B[i] - tmp;
    }
}

//UX = Y
void UX_Y(float** upper, float* Y, float* result_tab, int& size)
{
    float tmp;
    for (int i = size - 1; i >= 0; i--)
    {
        tmp = 0;
        for (int j = i + 1; j < size; j++)
        {
            tmp += upper[i][j] * result_tab[j];
        }
        result_tab[i] = (Y[i] - tmp) / upper[i][i];
    }
}

int main()
{

    int n = 0;

    //matrixes
    float** A;
    float** U;
    float** lower;
    float** upper;

    //tables
    float* B;
    float* X;
    float* Y;
    float* result_tab;

    cout << "Podaj wielkosc macierzy A: ";
    cin >> n;

    A = get_matrix(n);
    B = get_tab(n);
    U = get_matrix(n + 1);
    lower = get_matrix(n);
    upper = get_matrix(n);


    cout << "Podaj wspolczynniki macierzy\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cin >> A[i][j];
        cin >> B[i];
    }

    //GAUSS
    cout << "\nGauss\n";
    //matrix of coefficients
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++) {
            if (j != n)
            {
                U[i][j] = A[i][j];
            }
            else
            {
                U[i][j] = B[i];
            }
        }
    }
    cout << "===Before elimination===\n";
    show_coefficient_matrix(U, n);
    cout << "\n======\n";
//    partial_pivoting(matrix, n);
    forward_elim(U, n);
    cout << "===After elimination===\n";
    show_coefficient_matrix(U, n);
    result_tab = get_tab(n);
    backward_substitution(U, result_tab, n);
    for (int i = 0; i < n; i++)
        cout << "x" << i + 1 << ": " << result_tab[i] << endl;
    //END OF GAUSS SOLUTION

    //LU
    cout << "\nLU decomposition\n";
    Y = get_tab(n);
    X = get_tab(n);
    
    //put 1.0 on lower_triangular_matrix diagonal and reset result_tab[]
    for (int i = 0; i < n; i++)
    {
        lower[i][i] = 1.0;
        result_tab[i] = 0.0;
    }

    //LU Decomposition
    for (int i = 0; i < n; i++)
    {
        upper_triangular_matrix(A, upper, lower, n, i);
        if (i < (n - 1))
        {
            lower_triangular_matrix(A, upper, lower, n, i);
        }
    }
    cout << "===Upper triangular matrix===\n";
    show_matrix(upper, n);
    cout << "\n===Upper triangular matrix===\n";
    show_matrix(lower, n);

    //Solving L*U*X=B

    //ad.1 L*Y = B
    LY_B(lower, B, Y, n);
    //ad.2 U*X = Y
    UX_Y(upper, Y, result_tab, n);

    for (int i = 0; i < n; i++)
        cout << "x" << i + 1 << ": " << result_tab[i] << endl;

    //END OF LU SOLUTION

    return 0;

}