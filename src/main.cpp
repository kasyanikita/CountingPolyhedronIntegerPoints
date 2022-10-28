#include "global_defs.h"

int main()
{
    GroupIP::Timer t;
    std::vector<int> v = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    t.start();
    std::sort(v.begin(), v.end());
    t.finish() << "\n";
}