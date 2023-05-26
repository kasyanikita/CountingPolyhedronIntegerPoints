#include "include/SNF.h"
#include <iostream>

using namespace MSVP;

void test_my_DF(size_t retries = 100, size_t n = 100, unsigned int interval = 100)
{
    default_random_engine eng;
    int error_number = 0;
    for (size_t retr_i = 0; retr_i < retries; retr_i++)
    {
        int_mat A = generate_random_mat(n, interval, eng);
        int_mat P, S, Q;
        tie(P, S, Q) = my_DF(A);

        if (!is_diagonal(S))
        {
            cout << "S is not diagonal: " << endl;
            //cout << S << endl;
        }

        if (P * A * Q != S)
        {
            error_number++;
            cout << "DF computation contains errors. P A Q is: " << endl;
            //cout << P * A * Q << endl;

            //cout << "S is: " << endl;
            //cout << S << endl;
        }
    }
    cout << error_number << endl;
}

int main()
{
    test_my_DF(1000, 10, 10);
    return 0;
}