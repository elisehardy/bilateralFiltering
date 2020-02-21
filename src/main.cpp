

#include <bilateral/Image8u1.hpp>
#include <bilateral/Bilateral.hpp>


using namespace bilateral;


int main(int argc, char **argv) {
    Image8u1 src = Image8u1(argv[1]), dst;
    Bilateral::naive(src, dst, 12., 16.);
    dst.save("output.png");
    
    return 0;
}
