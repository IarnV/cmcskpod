#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>
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

int main(void) {
    unsigned seed = time(0), dim_num = 3;
    long double value = 0.0, v = 1.0;
    int i = 0;
    std::vector<std::pair<double,double> > borders = {std::make_pair(-3.56, -0.49), std::make_pair(-2.15, 0.72), std::make_pair(-6.13, -4.84)};
    for (std::vector<std::pair<double,double> >::iterator iter = borders.begin(); iter != borders.end(); iter++) {
        v *= std::abs((*iter).first - (*iter).second);
    }
    gener rng(1664525, 1013904223, seed, 1);
    unsigned iter_num = 100000000, timer = std::time(nullptr);
    for (i=0; i<iter_num; i++) {
        std::vector<double> x(dim_num);
        for (int j=0; j < dim_num; j++) {
            x[j] = borders[j].first + (borders[j].second - borders[j].first) * rng.next(0);

        }
        double y = f(x);
        value += y;
    }
    timer = std::time(nullptr) - timer;
    std::cout << "Время выполнения программы = " << timer << 's' << std::endl;
    double res = value / iter_num * v;
    std::cout << "Результат вычисления интеграла = " << res << std::endl;
    return 0;
}



    
