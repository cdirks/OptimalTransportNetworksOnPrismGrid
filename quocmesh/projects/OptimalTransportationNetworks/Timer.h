#ifndef Timer_h
#define Timer_h

#include<omp.h>

class Timer {
public:
    Timer(void) : t_begin() {};

    void reset() {
        t_begin = omp_get_wtime();
    }

    void end() {
        t_end = omp_get_wtime();
    }

    double elapsed() const {
        return t_end - t_begin;
    }

private:
    double t_begin;
    double t_end;
};

#endif