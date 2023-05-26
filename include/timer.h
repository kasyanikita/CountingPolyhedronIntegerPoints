#ifndef _TIMER_H_
#define _TIMER_H_
#include <chrono>
#include <sstream>

#include "boost_tee_logging.h"

namespace GroupIP
{
    using namespace std;

    class Timer
    {
    private:
        mutable chrono::steady_clock::time_point _start;

    public:
        void start() const
        {
            _start = chrono::steady_clock::now();
        }

        string finish_str(string suff = "'s") const
        {
            stringstream res;

            auto end = std::chrono::steady_clock::now();

            std::chrono::duration<double> elapsed_seconds = end - _start;

            res << elapsed_seconds.count() << suff;

            return res.str();
        }

        template <typename streamT = logger::stream>
        streamT &finish(streamT &out_stream = default_log.out(), string suff = "'s") const
        {
            out_stream << finish_str(suff);
            return out_stream;
        }
    };

    Timer global_timer;
} // namespace GrouIP

#endif