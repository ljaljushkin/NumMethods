#ifndef TIMER_H
#define TIMER_H


/* ������������������ ������ ��� Windows- � �������� *nix-������.
�������������:
... [��������, ������� �� ���� �����������: ������ �� ����� � �.�.]
Timer timer;
timer.start();
... [��������, ������� ���� �����������, - �������������� ����� ���������]
timer.stop();
double timeInSeconds = timer.getElapsed();
... [timeInSeconds - ���������� ����� � ��������, ��� ���� ������� � ����]
*/

#ifdef _WIN32

#include <windows.h>

class Timer
{
public:

    Timer() { reset(); }

    void reset() { total = 0; }

    void start() { QueryPerformanceCounter(&last_start); }

    void stop()
    {
        LARGE_INTEGER now;
        QueryPerformanceCounter(&now);
        total += now.QuadPart - last_start.QuadPart;
    }

    double getElapsed()
    {
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        return total / (double) freq.QuadPart;
    }

private:

    LARGE_INTEGER last_start;
    long long total;
};

#else

#include <time.h>

class Timer
{
public:

    Timer() { reset(); }
        
    void reset() { total = 0; }
        
    void start() { clock_gettime(CLOCK_MONOTONIC, &last_start); }
        
    void stop() {
        timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        total += (now.tv_sec - last_start.tv_sec) * 1000000000LL
            + now.tv_nsec - last_start.tv_nsec;
    }
        
    double getElapsed() { return total / 1e9; }

private:

    timespec last_start;
    long long total;
};

#endif


#endif
