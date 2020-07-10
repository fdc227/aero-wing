#include "main.hpp"
#include <fstream>
#include "mkl_lapacke.h"
#include <vector>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#define N 63
#define NRHS 1
#define LDA N
#define LDB NRHS
#define STORAGE_SIZE 64*10000
#define VEC_LEN 64

using namespace std;
using namespace boost::numeric::odeint;

vector<double> A(LDA * N);
vector<double> rhs(N);
vector<double> var_list(2*N);
double storage[STORAGE_SIZE];
int location = 0;
int total_array_num = 0;
int vector_length = VEC_LEN;


void init_condition(vector<double>& var_list)
{
    for (int i = 0; i < 21; i++)
    {
        if (i < 10)
        {
            var_list[3 * i] = 1.0 - 0.1 * i;
            var_list[3 * i + 1] = -0.1;
            var_list[3 * i + 2] = 1.0 - 0.1 * i;
        }
        else if (i == 10)
        {
            var_list[30] = 0;
            var_list[31] = 0;
            var_list[32] = 0;
        }
        else
        {
            var_list[3 * i] = 0.1 * (i - 10);
            var_list[3 * i + 1] = 0.1;
            var_list[3 * i + 2] = 0.1 * (i - 10);
        }
    }
}

void print_array(vector<double>& array, int length)
{
    for(int i=0;i<length;i++)
    {
        if(i!=length-1)
        cout<<array[i]<<' ';
        else
        {
            cout<<array[i]<<endl;
        }
        
    }
}

void array_from_vector(double* array_ptr, const vector<double>& V, const double t)
{
    array_ptr[0] = t;
    for (int i = 0; i < V.size(); i++)
    {
        array_ptr[i+1] = V[i];
    }
}

void ptr_array_to_console(const vector<double>& x, const double t)
{
    cout << t<< ' ';
    for (int i = 0; i < 63; i++)
    {
        cout << x[i] << ' ';
    }
    cout << endl;
}

void write_array_to_file(double* array, int total_array_num, int array_size, ofstream& file)
{
    for(int i=0;i<total_array_num;i++)
    {
        for (int j = 0; j < array_size; j++)
        {
            if (j != array_size - 1)
            {
                file<<array[array_size*i+ j]<<' ';
            }
            else
            {
                file<<array[array_size * i + j]<<endl;
            }
        }
    }
    for (int i = 0; i < total_array_num * array_size; i++)
    {
        array[i] = 0;
    }
    location = 0;
    total_array_num = 0;
}

void output_process(const vector<double>& x, const double t)
{
    total_array_num += 1;
    array_from_vector(&storage[location], x, t);
    if (location + vector_length >= STORAGE_SIZE)
    {
        ofstream outfile;
        outfile.open("ODE_RESULTS.dat", ios_base::app);
        write_array_to_file(storage, total_array_num, vector_length, outfile);
    }
    else
    {
        location += vector_length;
    }
    
}

void ODE_dydt(const vector<double>& var_list, vector<double>& rhs, double t)
{
    MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
    MKL_INT ipiv[N];

    A_func(var_list, A);
    rhs_func(var_list, rhs);
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, &*A.begin(), lda, ipiv, &*rhs.begin(), ldb);

}

int main(void)
{

    init_condition(var_list);

    integrate(ODE_dydt, var_list, 0.0, 10.0, 0.1, ptr_array_to_console);
   /* ofstream outfile;
    outfile.open("ODE_RESULTS.dat", ios_base::app);
    write_array_to_file(storage, total_array_num, vector_length, outfile);*/

    //A_func(var_list, A);
    //rhs_func(var_list, rhs);

    //info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, & *A.begin(), lda, ipiv, & *rhs.begin(), ldb);
    //print_array(rhs, 63);

    //ofstream outfile;
    //outfile.open("data.dat");*/

    //rhs_func(var_list, rhs);
    //print_array(rhs, 63);

    //save_array(rhs, 63, outfile);

    //outfile.close();

    return 0;

}