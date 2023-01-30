#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <mpi.h>
#include <cstdlib>

enum
{
    C_INT = 4294967295,
};

class gener {
private:
    std::vector<unsigned> seed;
    unsigned a;
    unsigned c;
public:
    gener(unsigned _a, unsigned _c, unsigned _seed, int s): a(_a), c(_c) { 
        seed.resize(s);
        for (std::vector<unsigned>::iterator it = seed.begin(); it != seed.end(); it++) {
            *it = (std::rand() + _seed) % C_INT;
        }
    };
    double next(int seed_num) {
        double ans = ((seed[seed_num] * a + c) % C_INT) / 4294967295.0;
        seed[seed_num] = (seed[seed_num] * a + c) % C_INT;
        return ans;
    };
};

double f(std::vector<double> x) {
    double ans = 0.0;
    for (std::vector<double>::iterator it = x.begin(); it != x.end(); it++) {
        ans += *it * *it;
    }
    return cos(-ans) + 7 * sin(ans);
}

int main(int argc, char *argv[]) {
    long double v = 1.0, sum_val = 0.0, res = 0.0;
    std::vector<std::pair<double,double> > borders = {std::make_pair(-3.56, -0.49), std::make_pair(-2.15, 0.72), std::make_pair(-6.13, -4.84)};
    int my_id, proc_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    unsigned seed = std::rand(), iter_num = 100000000, dim_num = 3;
    long double value = 0.0;
    double starttime, endtime;
    starttime = MPI_Wtime();
    gener rng(1664525, 1013904223, seed + 17 * my_id, 1);
    std::vector<double> x(dim_num);
    double vals = 0.0;
    for (int i = 0; i < iter_num / proc_num; i++) {
        for (int k = 0; k < dim_num; k++) {
            x[k] = borders[k].first + (borders[k].second - borders[k].first) * rng.next(0);
        }
        value += f(x);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&value, &vals, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_id == 0) {
        for (std::vector<std::pair<double,double> >::iterator iter = borders.begin(); iter != borders.end(); iter++) {
            v *= std::abs((*iter).first - (*iter).second);
        }
        res = vals / iter_num * v;
        endtime = MPI_Wtime();
        std::cout << "Время работы = " << endtime-starttime << "s" << std::endl;
        std::cout << "Результат вычисления интеграла = " << res << std::endl;
    }
    MPI_Finalize();
    return 0;
}
